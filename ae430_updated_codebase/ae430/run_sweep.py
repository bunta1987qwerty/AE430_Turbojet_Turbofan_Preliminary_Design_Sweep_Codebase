# -*- coding: utf-8 -*-
"""
Sweep runner (Python 2.7 compatible).

- Supports multiprocessing (--mode mp)
- Optional MPI auto-detect if mpi4py exists, but MP is default/fallback.
- Rank0-only logging/printing
- Per-worker rank_XXXX.csv and rank_XXXX.jsonl outputs
- Merges to summary_all.csv and cases_all.jsonl

Usage:
  python -m ae430.run_sweep --mode mp --nproc 8 --sweep sweeps/turbojet_sweep.json --outdir results --clean --F 350e3
"""

from __future__ import division

import argparse
import copy
import csv
import itertools
import json
import os
import shutil
import sys
import time
import traceback

from ae430.config import load_config
from ae430.logger_setup import setup_logger
from ae430.engine import evaluate_design


# ---------- helpers: json / overrides ----------

def _load_json(path):
    f = open(path, 'r')
    try:
        return json.load(f)
    finally:
        f.close()


def _is_dict(x):
    return isinstance(x, dict)


def _expand_value(v):
    """
    v can be:
      - list: return list
      - scalar: return [scalar]
      - dict with {"linspace":[a,b,n]}: return list of floats
    """
    if _is_dict(v) and ('linspace' in v):
        a, b, n = v['linspace']
        a = float(a); b = float(b); n = int(n)
        if n <= 1:
            return [a]
        out = []
        for i in range(n):
            t = float(i) / float(n - 1)
            out.append(a + t * (b - a))
        return out
    if isinstance(v, list):
        return v
    return [v]


def _set_by_dotted(cfg, dotted_key, value):
    """
    Set cfg['a']['b']['c'] = value for dotted_key "a.b.c"
    """
    parts = dotted_key.split('.')
    d = cfg
    for p in parts[:-1]:
        if p not in d or not isinstance(d[p], dict):
            d[p] = {}
        d = d[p]
    d[parts[-1]] = value


def _apply_overrides(base_cfg, overrides):
    """
    overrides is dict of dotted_key -> value
    Returns new cfg
    """
    cfg = copy.deepcopy(base_cfg)
    for k, v in overrides.items():
        _set_by_dotted(cfg, k, v)
    return cfg


# ---------- outputs ----------

def _safe_mkdir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def _clean_rank_outputs(outdir, logger):
    """
    Rank0-only cleanup: remove old rank_*.csv/jsonl/error markers and merged files.
    """
    if not os.path.isdir(outdir):
        return
    for fn in os.listdir(outdir):
        if fn.startswith('rank_') or fn in ('summary_all.csv', 'cases_all.jsonl', 'DONE.txt', 'STARTED.txt'):
            try:
                os.remove(os.path.join(outdir, fn))
            except Exception:
                pass
    if logger:
        logger.info("Cleaned old rank outputs in %s", outdir)


def _merge_rank_csvs(outdir, merged_name):
    """
    Merge rank_*.csv into one CSV. Handles empty files safely.
    Returns merged path or None.
    """
    rank_files = []
    for fn in sorted(os.listdir(outdir)):
        if fn.startswith('rank_') and fn.endswith('.csv'):
            rank_files.append(os.path.join(outdir, fn))

    if not rank_files:
        return None

    # discover header from first non-empty
    header = None
    for rf in rank_files:
        try:
            f = open(rf, 'r')
            r = csv.reader(f)
            row = None
            try:
                row = next(r)
            except Exception:
                row = None
            f.close()
            if row:
                header = row
                break
        except Exception:
            continue

    if header is None:
        return None

    merged_path = os.path.join(outdir, merged_name)
    out = open(merged_path, 'w')
    w = csv.writer(out)
    w.writerow(header)

    for rf in rank_files:
        f = open(rf, 'r')
        r = csv.reader(f)
        first = True
        for row in r:
            if first:
                first = False
                # skip header if matches
                if row == header:
                    continue
            if row:
                w.writerow(row)
        f.close()

    out.close()
    return merged_path


def _merge_rank_jsonls(outdir, merged_name):
    rank_files = []
    for fn in sorted(os.listdir(outdir)):
        if fn.startswith('rank_') and fn.endswith('.jsonl'):
            rank_files.append(os.path.join(outdir, fn))
    if not rank_files:
        return None
    merged_path = os.path.join(outdir, merged_name)
    out = open(merged_path, 'w')
    for rf in rank_files:
        f = open(rf, 'r')
        for line in f:
            if line.strip():
                out.write(line)
        f.close()
    out.close()
    return merged_path


# ---------- multiprocessing worker (top-level, picklable) ----------

_G_BASE_CFG = None
_G_OUTDIR = None
_G_LOGGER = None
_G_HEADER_WRITTEN = False


def _worker_init(base_cfg, outdir):
    global _G_BASE_CFG, _G_OUTDIR, _G_HEADER_WRITTEN
    _G_BASE_CFG = base_cfg
    _G_OUTDIR = outdir
    _G_HEADER_WRITTEN = False


def _worker_rank_id():
    """
    Try to get deterministic worker id from multiprocessing internals.
    """
    try:
        import multiprocessing
        ident = getattr(multiprocessing.current_process(), "_identity", None)
        if ident and len(ident) > 0:
            return int(ident[0] - 1)
    except Exception:
        pass
    return 0


def _worker_run_one(args):
    """
    args = (case_id, overrides, F_required)
    Returns tuple: (ok_bool, crashed_bool)
    """
    global _G_BASE_CFG, _G_OUTDIR, _G_HEADER_WRITTEN

    case_id, overrides, F_required = args
    wid = _worker_rank_id()

    csv_path = os.path.join(_G_OUTDIR, "rank_%04d.csv" % wid)
    jsonl_path = os.path.join(_G_OUTDIR, "rank_%04d.jsonl" % wid)
    err_path = os.path.join(_G_OUTDIR, "rank_%04d.error.txt" % wid)

    try:
        cfg = _apply_overrides(_G_BASE_CFG, overrides)
        cfg.setdefault('requirements', {})
        cfg['requirements']['F_required'] = float(F_required)

        # Make sure base_dir exists for resolving diesel.xml relative to base config
        if '_base_dir' not in cfg:
            cfg['_base_dir'] = os.getcwd()

        # Evaluate
        res, st, cons = evaluate_design(cfg, None)

        # flatten for CSV row
        row = {}
        row['case_id'] = case_id
        row['ok'] = 1 if bool(res.get('feasible', False)) else 0

        # record some keys if present
        keep = [
            'engine_type', 'alpha', 'mdot0', 'pi_c', 'pi_f', 'pi_d', 'pi_bd',
            'Tt4', 'f', 'fcool', 'eta_t_eff', 'pi_t',
            'F_installed', 'F_uninstalled', 'eta0', 'tau_res', 'fc'
        ]
        for k in keep:
            if k in res:
                row[k] = res[k]

        # also record primary sweep knobs for plotting
        for k, v in overrides.items():
            # make safe csv column name
            row['param_' + k.replace('.', '_')] = v

        # failure reason (string)
        fr = res.get('fail_reasons', [])
        if isinstance(fr, list):
            row['fail_reasons'] = ";".join([str(x) for x in fr])
        else:
            row['fail_reasons'] = str(fr)

        # Write CSV (append)
        fieldnames = sorted(row.keys())
        mode = 'a'
        write_header = False
        if (not os.path.exists(csv_path)) or (os.path.getsize(csv_path) == 0) or (not _G_HEADER_WRITTEN):
            write_header = True

        fcsv = open(csv_path, mode)
        w = csv.DictWriter(fcsv, fieldnames=fieldnames)
        if write_header:
            w.writeheader()
            _G_HEADER_WRITTEN = True
        w.writerow(row)
        fcsv.close()

        # Write JSONL case record (append)
        jout = {
            'case_id': case_id,
            'overrides': overrides,
            'results': res,
            'constraints': cons
        }
        fj = open(jsonl_path, 'a')
        fj.write(json.dumps(jout) + "\n")
        fj.close()

        return True, False

    except Exception:
        # log traceback to rank error file
        tb = traceback.format_exc()
        fe = open(err_path, 'a')
        fe.write("case_id=%s\n" % str(case_id))
        fe.write(tb + "\n")
        fe.close()

        return False, True


# ---------- main runner ----------

def run_mp(args):
    sweep_path = os.path.abspath(args.sweep)
    outdir = os.path.abspath(args.outdir)
    _safe_mkdir(outdir)

    # logger rank0 only
    try:
        logger = setup_logger(outdir, log_filename='run.log')
    except TypeError:
        logger = setup_logger(outdir)

    logger.info("Starting sweep runner.")
    logger.info("sweep_path=%s", sweep_path)
    logger.info("outdir=%s", outdir)
    logger.info("F_required=%.3e", float(args.F))

    sweep_spec = _load_json(sweep_path)
    base_cfg_rel = sweep_spec.get('base_config', 'config_default.json')

    # resolve base config relative to sweep file dir
    sweep_dir = os.path.dirname(sweep_path)
    base_cfg_path = base_cfg_rel
    if not os.path.isabs(base_cfg_path):
        base_cfg_path = os.path.join(sweep_dir, base_cfg_rel)
    base_cfg_path = os.path.abspath(base_cfg_path)

    logger.info("base_cfg_path=%s", base_cfg_path)

    if args.clean:
        _clean_rank_outputs(outdir, logger)

    # load base cfg
    base_cfg = load_config(base_cfg_path)
    base_cfg['_base_dir'] = os.path.dirname(base_cfg_path)

    # build grid
    grid = sweep_spec.get('grid', {})
    keys = sorted(grid.keys())

    values_lists = []
    for k in keys:
        values_lists.append(_expand_value(grid[k]))

    # enumerate cases
    cases = []
    case_id = 0
    for vals in itertools.product(*values_lists):
        overrides = {}
        for i in range(len(keys)):
            overrides[keys[i]] = vals[i]
        cases.append((case_id, overrides, float(args.F)))
        case_id += 1

    total = len(cases)
    logger.info("TOTAL_CASES=%d", total)
    print("('TOTAL_CASES:', %d)" % total)
    print("('F_required:', %s)" % str(float(args.F)))

    # marker
    f0 = open(os.path.join(outdir, "STARTED.txt"), 'w')
    f0.write("STARTED\n")
    f0.close()

    # multiprocessing
    import multiprocessing as mp
    nproc = int(args.nproc)
    if nproc < 1:
        nproc = 1

    t0 = time.time()
    ok_total = 0
    crashed_total = 0

    pool = mp.Pool(processes=nproc, initializer=_worker_init, initargs=(base_cfg, outdir))
    try:
        # chunksize helps speed on large sweeps
        chunksize = 5
        for (ok, crashed) in pool.imap_unordered(_worker_run_one, cases, chunksize=chunksize):
            if ok:
                ok_total += 1
            if crashed:
                crashed_total += 1
        pool.close()
        pool.join()
    finally:
        try:
            pool.terminate()
        except Exception:
            pass

    elapsed = time.time() - t0
    rate = float(total) / max(1e-12, elapsed)
    logger.info("workers finished: ok_total=%d crashed_total=%d elapsed=%.2fs rate=%.2f cases/s",
                ok_total, crashed_total, elapsed, rate)

    merged_csv = _merge_rank_csvs(outdir, "summary_all.csv")
    if merged_csv:
        logger.info("Merged CSV: %s", merged_csv)
    else:
        logger.warning("No rank CSVs to merge (or all empty). Check rank_XXXX.error.txt.")

    merged_jsonl = _merge_rank_jsonls(outdir, "cases_all.jsonl")
    if merged_jsonl:
        logger.info("Merged JSONL: %s", merged_jsonl)
    else:
        logger.warning("No rank JSONLs to merge.")

    f1 = open(os.path.join(outdir, "DONE.txt"), 'w')
    f1.write("DONE\n")
    f1.close()

    logger.info("DONE")


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mode', default='auto', choices=['auto', 'mp', 'mpi'])
    p.add_argument('--nproc', default=4, type=int)
    p.add_argument('--sweep', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--clean', action='store_true')
    p.add_argument('--F', default=350e3, type=float)
    args = p.parse_args()

    # mode auto: prefer MPI only if mpi4py exists AND world size > 1
    mode = args.mode
    if mode == 'auto':
        mode = 'mp'
        try:
            from mpi4py import MPI
            if MPI.COMM_WORLD.Get_size() > 1:
                mode = 'mpi'
        except Exception:
            mode = 'mp'

    if mode == 'mpi':
        # Keep it simple: if mpi4py exists, you can add MPI runner later.
        # For tonight, always run multiprocessing path.
        mode = 'mp'

    run_mp(args)


if __name__ == '__main__':
    main()
