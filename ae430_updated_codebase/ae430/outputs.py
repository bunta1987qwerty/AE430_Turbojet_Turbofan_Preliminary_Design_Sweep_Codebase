# -*- coding: utf-8 -*-
from __future__ import division

import csv
import json
import os
import glob

try:
    text_types = (str, unicode)
except NameError:
    text_types = (str,)

def _safe_float(x):
    try:
        return float(x)
    except Exception:
        return x

def write_jsonl_line(path, obj):
    with open(path, 'a') as f:
        f.write(json.dumps(obj))
        f.write("\n")

def write_rank_csv_header_if_needed(path, fieldnames):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return
    with open(path, 'w') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()

def append_rank_csv_row(path, fieldnames, rowdict):
    out = {}
    for k in fieldnames:
        out[k] = rowdict.get(k, '')
    with open(path, 'a') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writerow(out)

def merge_rank_csvs(outdir, merged_name='summary_all.csv'):
    paths = sorted(glob.glob(os.path.join(outdir, "rank_*.csv")))
    # ignore empty
    paths = [p for p in paths if os.path.getsize(p) > 0]
    if not paths:
        return None

    # find a header
    header = None
    for p in paths:
        with open(p, 'r') as f:
            r = csv.DictReader(f)
            if r.fieldnames:
                header = list(r.fieldnames)
                break
    if not header:
        return None

    merged_path = os.path.join(outdir, merged_name)
    with open(merged_path, 'w') as out_f:
        w = csv.DictWriter(out_f, fieldnames=header)
        w.writeheader()
        for p in paths:
            with open(p, 'r') as in_f:
                r = csv.DictReader(in_f)
                if not r.fieldnames:
                    continue
                for row in r:
                    w.writerow(row)
    return merged_path

def merge_rank_jsonls(outdir, merged_name='cases_all.jsonl'):
    paths = sorted(glob.glob(os.path.join(outdir, "rank_*.jsonl")))
    paths = [p for p in paths if os.path.getsize(p) > 0]
    if not paths:
        return None
    merged_path = os.path.join(outdir, merged_name)
    with open(merged_path, 'w') as out_f:
        for p in paths:
            with open(p, 'r') as in_f:
                for line in in_f:
                    out_f.write(line)
    return merged_path
