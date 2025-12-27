#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run the AE430 final project pipeline.

Usage (inside your AE430 conda env):
    python run_pipeline.py --config config_default.json

Optional sweep:
    python run_pipeline.py --config config_default.json --sweep_pi_c 5 60 30

Outputs:
    outputs/<run_name>/
        run.log
        results.json
        results.csv
        stations.csv
        constraints.csv
        mass_flows.png  (optional)
        sweep_pi_c.csv  (optional)
        sweep_pi_c.png  (optional)
"""
from __future__ import division
import argparse
import copy
import os
import time

from ae430.config import load_config, validate_config
from ae430.logger_setup import setup_logger
from ae430.engine import evaluate_design
from ae430.outputs import write_json, write_results_csv, write_stations_csv, write_constraints_csv, make_basic_plots

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True, help='Path to JSON config file')
    p.add_argument('--out_dir', default=None, help='Output directory (default: from config or outputs/<timestamp>)')
    p.add_argument('--debug', action='store_true', help='Enable DEBUG logging')
    p.add_argument('--no_plots', action='store_true', help='Disable plotting')
    p.add_argument('--sweep_pi_c', nargs=3, metavar=('PI_C_MIN','PI_C_MAX','N'),
                   help='Sweep compressor pressure ratio pi_c over a range (inclusive)')
    return p.parse_args()

def main():
    args = parse_args()
    cfg = load_config(args.config)

    missing = validate_config(cfg)
    if missing:
        raise SystemExit("Missing required config fields: %s" % ', '.join(missing))

    out_base = cfg.get('outputs', {}).get('out_dir', 'outputs')
    if args.out_dir:
        out_dir = args.out_dir
    else:
        stamp = time.strftime('%Y%m%d_%H%M%S')
        out_dir = os.path.join(out_base, 'run_' + stamp)

    log_level = 10 if args.debug else 20
    logger = setup_logger(out_dir, log_filename=cfg.get('outputs', {}).get('log_file', 'run.log'),
                          level=log_level)

    logger.info("Config: %s", cfg.get('_config_path'))

    # Single design evaluation
    results, stations, constraints = evaluate_design(cfg, logger)

    write_json(results, os.path.join(out_dir, 'results.json'))
    write_results_csv(results, os.path.join(out_dir, 'results.csv'))
    write_stations_csv(stations, os.path.join(out_dir, 'stations.csv'))
    write_constraints_csv(constraints, os.path.join(out_dir, 'constraints.csv'))

    if (not args.no_plots) and bool(cfg.get('outputs', {}).get('make_plots', True)):
        make_basic_plots(results, out_dir)

    # Optional sweep over pi_c
    if args.sweep_pi_c:
        pi_c_min = float(args.sweep_pi_c[0])
        pi_c_max = float(args.sweep_pi_c[1])
        N = int(args.sweep_pi_c[2])
        if N < 2:
            raise SystemExit("N must be >= 2")

        logger.info("Running pi_c sweep: %.2f -> %.2f (%d points)", pi_c_min, pi_c_max, N)

        sweep_rows = []
        for k in range(N):
            pi_c = pi_c_min + (pi_c_max - pi_c_min) * (k / float(N-1))
            cfg_k = copy.deepcopy(cfg)
            cfg_k['design']['pi_c'] = pi_c
            try:
                res_k, _, cons_k = evaluate_design(cfg_k, logger)
            except Exception as e:
                logger.warning("pi_c=%.2f failed: %s", pi_c, str(e))
                continue
            sweep_rows.append({
                'pi_c': pi_c,
                'fc': res_k.get('fc', None),
                'eta0': res_k.get('eta0', None),
                'mdot0': res_k.get('mdot0', None),
                'Tt4': res_k.get('Tt4', None),
                'XCO': res_k.get('XCO', None),
            })

        # Write sweep CSV
        import csv
        sweep_path = os.path.join(out_dir, 'sweep_pi_c.csv')
        with open(sweep_path, 'w') as f:
            w = csv.DictWriter(f, fieldnames=['pi_c','fc','eta0','mdot0','Tt4','XCO'])
            w.writeheader()
            for r in sweep_rows:
                w.writerow(r)
        logger.info("Wrote sweep results to %s", sweep_path)

        # Plot sweep if possible
        if (not args.no_plots) and bool(cfg.get('outputs', {}).get('make_plots', True)):
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                pis = [r['pi_c'] for r in sweep_rows]
                fcs = [r['fc'] for r in sweep_rows]
                plt.figure()
                plt.plot(pis, fcs, marker='o')
                plt.xlabel('Compressor pressure ratio, pi_c')
                plt.ylabel('Cost function f_c')
                plt.title('f_c vs pi_c sweep')
                plt.tight_layout()
                plt.savefig(os.path.join(out_dir, 'sweep_pi_c.png'))
                plt.close()
            except Exception as e:
                logger.warning("Could not plot sweep: %s", str(e))

    logger.info("Done. Outputs in %s", out_dir)

if __name__ == '__main__':
    main()
