# -*- coding: utf-8 -*-
from __future__ import division

import argparse
import csv
import math

def _to_float(x):
    try:
        return float(x)
    except Exception:
        return None

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", required=True, help="Path to summary_all.csv")
    p.add_argument("--x", required=True, help="Column name for x axis")
    p.add_argument("--y", required=True, help="Column name for y axis")
    p.add_argument("--z", required=True, help="Column name for z axis")
    p.add_argument("--c", required=True, help="Column name for color")
    p.add_argument("--filter-feasible", action="store_true", help="Keep only feasible==1")
    args = p.parse_args()

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    xs, ys, zs, cs = [], [], [], []

    with open(args.csv, "r") as f:
        r = csv.DictReader(f)
        for row in r:
            if args.filter_feasible:
                feas = row.get("feasible", "0")
                if feas not in ("1", "True", "true"):
                    continue

            x = _to_float(row.get(args.x))
            y = _to_float(row.get(args.y))
            z = _to_float(row.get(args.z))
            c = _to_float(row.get(args.c))

            if x is None or y is None or z is None or c is None:
                continue
            if math.isnan(x) or math.isnan(y) or math.isnan(z) or math.isnan(c):
                continue

            xs.append(x); ys.append(y); zs.append(z); cs.append(c)

    if len(xs) == 0:
        print("No points found. Check column names and feasibility filter.")
        return

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(xs, ys, zs, c=cs)
    ax.set_xlabel(args.x)
    ax.set_ylabel(args.y)
    ax.set_zlabel(args.z)
    cb = plt.colorbar(sc)
    cb.set_label(args.c)
    plt.title("AE430 sweep 3D")
    plt.show()

if __name__ == "__main__":
    main()
