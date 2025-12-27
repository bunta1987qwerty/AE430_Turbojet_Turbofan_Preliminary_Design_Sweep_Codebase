# -*- coding: utf-8 -*-
from __future__ import division
import math

def turbine_power_match(Tt4, pt4, w_req_per_kg_core, f, eta_t, gamma, cp):
    """
    Lumped turbine power match.

    Turbine sees mdot_gas = mdot_core*(1+f). If w_req is per kg core air, then
      w_turb_per_kg_gas = w_req_per_kg_core / (1+f)

    Turbine efficiency definition (total-to-total):
      eta_t = (Tt4 - Tt5) / (Tt4 - Tt5s)
    with
      (Tt5s/Tt4) = (pt5/pt4)^((gamma-1)/gamma)

    We choose Tt5 from work required, then infer pt5 from eta_t.
    """
    try:
        Tt4 = float(Tt4); pt4 = float(pt4)
        w_req = float(w_req_per_kg_core)
        f = float(f); eta_t = float(eta_t)
        gamma = float(gamma); cp = float(cp)
    except Exception:
        return {'ok': False, 'fail_reason': 'bad_inputs'}

    if (1.0 + f) <= 1e-9:
        return {'ok': False, 'fail_reason': 'bad_f'}

    w_per_kg_gas = w_req / (1.0 + f)
    dT = w_per_kg_gas / max(1e-12, cp)
    Tt5 = Tt4 - dT

    if Tt5 <= 50.0 or Tt5 >= Tt4:
        return {'ok': False, 'fail_reason': 'nonphysical_Tt5'}

    if eta_t <= 1e-6 or eta_t > 1.0:
        return {'ok': False, 'fail_reason': 'bad_eta_t'}

    # From eta definition:
    # Tt4 - Tt5s = (Tt4 - Tt5)/eta_t  -> Tt5s = Tt4 - (Tt4 - Tt5)/eta_t
    Tt5s = Tt4 - (Tt4 - Tt5) / eta_t
    if Tt5s <= 1e-9 or Tt5s >= Tt4:
        return {'ok': False, 'fail_reason': 'nonphysical_Tt5s'}

    # isentropic relation:
    # (pt5/pt4) = (Tt5s/Tt4)^(gamma/(gamma-1))
    expn = gamma / (gamma - 1.0)
    pi_t = (Tt5s / Tt4) ** expn
    pt5 = pt4 * pi_t

    if pi_t <= 0.0 or pi_t >= 1.0:
        # Turbine should drop pressure (pi_t < 1)
        return {'ok': False, 'fail_reason': 'bad_pi_t'}

    return {
        'ok': True,
        'Tt5': Tt5,
        'pt5': pt5,
        'pi_t': pi_t,
        'w_req_per_kg_core': w_req_per_kg_core,
        'w_per_kg_gas': w_per_kg_gas,
    }

# Backward-compat name some versions of engine.py used
def turbine_from_required_power(*args, **kwargs):
    return turbine_power_match(*args, **kwargs)
