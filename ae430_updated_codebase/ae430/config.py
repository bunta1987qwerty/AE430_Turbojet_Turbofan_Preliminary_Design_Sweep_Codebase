# -*- coding: utf-8 -*-
from __future__ import division

import copy
import json
import os

def load_config(path):
    with open(path, 'r') as f:
        cfg = json.load(f)
    return cfg

def deep_get(d, keypath, default=None):
    cur = d
    for k in keypath:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur

def deep_set(d, keypath, value):
    cur = d
    for k in keypath[:-1]:
        if k not in cur or not isinstance(cur[k], dict):
            cur[k] = {}
        cur = cur[k]
    cur[keypath[-1]] = value

def apply_overrides(cfg, overrides):
    """
    overrides: dict of dotted_key -> value
    Example: {"design.pi_c": 25.0, "losses.eta_c": 0.9}
    """
    out = copy.deepcopy(cfg)
    for dk, v in overrides.items():
        if isinstance(dk, unicode):
            dk = str(dk)
        parts = dk.split('.')
        deep_set(out, parts, v)
    return out

def validate_required(cfg):
    """
    Returns list of missing key strings (best-effort).
    """
    required = [
        ('flight', 'T0'),
        ('flight', 'p0'),
        ('flight', 'V0'),
        ('requirements', 'F_required'),
        ('design', 'engine_type'),
        ('design', 'pi_c'),
        ('losses', 'eta_m'),
        ('losses', 'eta_n'),
    ]
    missing = []
    for a, b in required:
        if cfg.get(a, {}).get(b, None) is None:
            missing.append("%s.%s" % (a, b))
    return missing

def resolve_path(base_dir, maybe_rel_path):
    if maybe_rel_path is None:
        return None
    if os.path.isabs(maybe_rel_path):
        return maybe_rel_path
    return os.path.abspath(os.path.join(base_dir, maybe_rel_path))
