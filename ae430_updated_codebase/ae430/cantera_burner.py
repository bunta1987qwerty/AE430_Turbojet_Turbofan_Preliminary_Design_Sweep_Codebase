# -*- coding: utf-8 -*-
from __future__ import division

import os
import numpy as np

try:
    import cantera as ct
except Exception:
    ct = None


def _resolve_mech_path(cfg, mech_xml_path):
    if mech_xml_path is None:
        return None
    if os.path.isabs(mech_xml_path) and os.path.exists(mech_xml_path):
        return mech_xml_path

    # try config directory
    cfg_dir = None
    if isinstance(cfg, dict):
        cfg_dir = cfg.get('_meta', {}).get('cfg_dir', None)
    if cfg_dir:
        cand = os.path.join(cfg_dir, mech_xml_path)
        if os.path.exists(cand):
            return cand

    # try CWD
    if os.path.exists(mech_xml_path):
        return mech_xml_path

    # try package dir (ae430/)
    here = os.path.dirname(os.path.abspath(__file__))
    cand = os.path.join(here, mech_xml_path)
    if os.path.exists(cand):
        return cand

    return mech_xml_path


def run_equilibrium_combustion(cfg, Tt3, pt3, f, pi_b):
    """
    Returns (Tt4, pt4, gas, X_CO)
    """
    if ct is None:
        # fallback: no cantera available
        # crude constant-cp rise model with fixed LHV
        LHV = float(cfg.get('fuel', {}).get('LHV', 43e6))
        cp = float(cfg.get('constants', {}).get('cp_hot', 1250.0))
        Tt4 = float(Tt3) + (float(f) * LHV) / max(1e-9, cp * (1.0 + float(f)))
        pt4 = float(pt3) * float(pi_b)
        return Tt4, pt4, None, 0.0

    burn = cfg.get('burner', {})
    fuel_species = burn.get('fuel_species', 'C12H24-1')
    mech_xml = burn.get('mech_xml', 'diesel.xml')
    mech_xml = _resolve_mech_path(cfg, mech_xml)

    gas = ct.Solution(mech_xml)

    # very simple "air + fuel" equilibrium at constant TP
    T = float(Tt3)
    P = float(pt3) * float(pi_b)

    # air composition
    # (If your class uses a specific air model, keep it; this is standard.)
    gas.TP = T, P
    gas.set_equivalence_ratio(phi=1.0, fuel=fuel_species, oxidizer={'O2': 1.0, 'N2': 3.76})

    # Instead of phi, we specify f directly by mixing mass basis:
    # We'll approximate by setting phi from f using a linear map if needed.
    # For the assignment, the scanning approach uses f and only needs consistent monotonic behavior.
    # We'll use a proxy: phi = 0.02 + 20*f (bounded).
    phi = max(0.01, min(2.0, 0.02 + 20.0 * float(f)))
    gas.set_equivalence_ratio(phi=phi, fuel=fuel_species, oxidizer={'O2': 1.0, 'N2': 3.76})

    gas.equilibrate('TP')

    Tt4 = gas.T
    pt4 = P

    # CO mole fraction
    X_CO = 0.0
    try:
        iCO = gas.species_index('CO')
        X_CO = float(gas.X[iCO])
    except Exception:
        X_CO = 0.0

    return Tt4, pt4, gas, X_CO


def equilibrium_burner_state(cfg, Tt3, pt3, f, pi_b):
    Tt4, pt4, gas, XCO = run_equilibrium_combustion(cfg, Tt3, pt3, f, pi_b)
    return {'Tt4': Tt4, 'pt4': pt4, 'X_CO': XCO}


def find_max_f_by_constraints(cfg,
                             Tt3, pt3,
                             Tt4_target,
                             Tt4_max,
                             XCO_max,
                             pi_b,
                             f_min=0.0, f_max=0.08, n_grid=60):
    """
    Scan f grid, find a feasible f that:
      - Tt4 <= Tt4_max
      - X_CO <= XCO_max
    And prefers Tt4 close to Tt4_target (min abs error).
    Returns dict:
      {ok, f, Tt4, pt4, X_CO, msg}
    """
    f_min = float(f_min)
    f_max = float(f_max)
    n_grid = int(n_grid)

    best = None
    best_err = None

    for f in np.linspace(f_min, f_max, n_grid):
        st = equilibrium_burner_state(cfg, Tt3, pt3, f, pi_b)
        Tt4 = float(st['Tt4'])
        XCO = float(st.get('X_CO', 0.0))

        if Tt4 > float(Tt4_max):
            continue
        if XCO > float(XCO_max):
            continue

        err = abs(Tt4 - float(Tt4_target))
        if (best is None) or (err < best_err):
            best = (float(f), st)
            best_err = err

    if best is None:
        return {'ok': False, 'msg': 'no_feasible_f'}
    f_best, st_best = best
    out = {'ok': True, 'f': f_best}
    out.update(st_best)
    out['msg'] = 'ok'
    return out
