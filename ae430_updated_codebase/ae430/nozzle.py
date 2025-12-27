# -*- coding: utf-8 -*-
"""
ae430/nozzle.py

Nozzle expansion model (simple, traceable, Python 2.7 compatible).

This file historically caused KeyError mismatches because some code expected
nozzle['F'] while others used nozzle['F_total'].

This version returns BOTH keys:
  - 'F_total' : momentum + pressure thrust (axial, incl. cosine loss)
  - 'F'       : alias of 'F_total' (for backward compatibility)

Assumptions:
  - Perfect expansion to ambient by default (p_exit = p_amb).
  - Isentropic relations with an efficiency eta_n applied to velocity.
"""
from __future__ import division, print_function

import math

from ae430.thermo import R_from_cp_gamma


def _nan():
    try:
        return float('nan')
    except Exception:
        return 0.0 / 0.0


def nozzle_expand_to_ambient(pt_in, Tt_in, p_amb, gamma, cp, eta_n, mdot,
                            half_angle_deg=7.0,
                            assume_perfect_expansion=True):
    """
    Expand from (pt_in, Tt_in) to p_amb (or compute no expansion if pt_in<=p_amb).

    Returns dict with:
        p_exit, T_exit, V_exit, A_exit,
        F_momentum, F_pressure, F_total, F
    """
    pt_in = float(pt_in)
    Tt_in = float(Tt_in)
    p_amb = float(p_amb)
    gamma = float(gamma)
    cp = float(cp)
    eta_n = float(eta_n)
    mdot = float(mdot)

    R = R_from_cp_gamma(cp, gamma)

    angle_factor = math.cos(float(half_angle_deg) * math.pi / 180.0)

    if not assume_perfect_expansion:
        raise NotImplementedError("Non-perfect-expansion nozzle not implemented in this simplified model.")

    # If pt_in <= p_amb, no useful expansion; set V ~ 0.
    if pt_in <= p_amb:
        p_exit = pt_in
        T_exit = Tt_in
        V = 0.0
    else:
        p_exit = p_amb
        # Isentropic relation for exit temperature
        T_exit_is = Tt_in * (p_exit / pt_in) ** ((gamma - 1.0) / gamma)
        # Apply nozzle efficiency to kinetic energy:
        dh_is = cp * (Tt_in - T_exit_is)
        if dh_is < 0.0:
            dh_is = 0.0
        V = math.sqrt(max(0.0, 2.0 * eta_n * dh_is))
        T_exit = T_exit_is

    rho_exit = p_exit / max(1.0e-12, (R * T_exit))
    A_exit = mdot / max(1.0e-12, (rho_exit * max(1.0e-9, V)))

    F_momentum = mdot * V * angle_factor
    F_pressure = (p_exit - p_amb) * A_exit * angle_factor
    F_total = F_momentum + F_pressure

    return {
        'p_exit': p_exit,
        'T_exit': T_exit,
        'V_exit': V,
        'A_exit': A_exit,
        'F_momentum': F_momentum,
        'F_pressure': F_pressure,
        'F_total': F_total,
        'F': F_total,  # alias for compatibility
        'angle_factor': angle_factor,
        'assume_perfect_expansion': bool(assume_perfect_expansion),
    }
