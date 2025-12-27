# -*- coding: utf-8 -*-
"""
Cooling model (Python 2.7 compatible).

Goal:
- Compute coolant fraction needed to satisfy Twg_max / Tc_max constraints.
- Convert uncooled turbine efficiency -> effective turbine efficiency penalty.

This file is intentionally conservative and "constraint-driven":
If cooling is enabled, we increase fcool until Twg <= Twg_max (and Tc <= Tc_max).

NOTES:
- We avoid relying on thermo.static_from_total_and_V in case thermo.py differs
  across teammates. We implement our own helper with the same physics.
"""

from __future__ import division
import math


def _R_from_cp_gamma(cp, gamma):
    return cp * (gamma - 1.0) / gamma


def _speed_of_sound(T, gamma, R):
    return math.sqrt(max(0.0, gamma * R * T))


def _static_from_total_and_V(Tt, pt, V, gamma, cp):
    """
    Constant-cp energy relation:
      cp*Tt = cp*T + V^2/2 => T = Tt - V^2/(2*cp)
    Then infer M from a(T), and p from isentropic relation.
    """
    R = _R_from_cp_gamma(cp, gamma)
    T = Tt - (V * V) / (2.0 * cp)
    if T <= 1.0:
        T = 1.0
    a = _speed_of_sound(T, gamma, R)
    M = V / a if a > 1e-12 else 0.0
    denom = (1.0 + 0.5 * (gamma - 1.0) * M * M) ** (gamma / (gamma - 1.0))
    p = pt / denom
    rho = p / (R * T)
    return T, p, rho, M


def solve_cooling_and_eta(cfg, Tt4, Tt3):
    """
    Parameters
    ----------
    cfg : dict
      expects:
        cfg['cooling'] and cfg['losses'] and cfg['requirements']
    Tt4 : float
      turbine inlet total temperature (station 4)
    Tt3 : float
      compressor exit total temperature (station 3) used as coolant source total temp

    Returns
    -------
    fcool : float
      coolant fraction (relative to core air)
    eta_t_eff : float
      effective turbine efficiency after cooling penalty
    ok : bool
    reason : str
    """
    cooling = cfg.get('cooling', {})
    req = cfg.get('requirements', {})
    losses = cfg.get('losses', {})

    enabled = bool(cooling.get('enabled', True))
    eta_t0 = float(losses.get('eta_t_uncooled', 0.90))

    if not enabled:
        return 0.0, eta_t0, True, "cooling_disabled"

    # Limits
    Twg_max = float(req.get('Twg_max', 1300.0))
    Tc_max = float(req.get('Tc_max', 1000.0))

    # Model inputs (simple 1D wall/film proxy)
    hg = float(cooling.get('hg', 7000.0))     # hot gas side h
    kw = float(cooling.get('kw', 50.0))       # wall conductivity
    tw = float(cooling.get('tw', 0.002))      # wall thickness
    Vc = float(cooling.get('Vc', 120.0))      # coolant speed proxy
    mu = float(cooling.get('mu', 4e-05))
    Pr = float(cooling.get('Pr', 0.7))

    # Use cold-side properties from constants, but if not present use typicals
    cst = cfg.get('constants', {})
    gamma_cold = float(cst.get('gamma_cold', 1.4))
    cp_cold = float(cst.get('cp_cold', 1000.0))

    # Coolant stagnation conditions (assume from compressor exit)
    # We don't have pt3 here; use a conservative assumption: coolant static derived from total only.
    # This is adequate for *constraint direction* (higher Vc -> lower static T).
    pt3 = 1.0  # dummy; pressure cancels for temperature checks in this simplified model

    # Compute coolant static inlet temperature proxy
    Tc_in, _, _, _ = _static_from_total_and_V(Tt3, pt3, Vc, gamma_cold, cp_cold)

    # Hot gas temperature near wall: assume Twg tends toward Tt4 (conservative)
    Tg = float(Tt4)

    # Simple conductive + convective wall model:
    # q" = hg*(Tg - Twg)
    # Twg - Twc = q"*tw/kw
    # q" = hc*(Twc - Tc) but we don't explicitly model hc; instead
    # use coolant fraction to reduce wall driving ΔT (proxy for film effectiveness).
    #
    # We'll model film effectiveness eps_f = 1 - exp(-K * fcool),
    # where eps_f increases with fcool. Then "effective gas temp" seen by wall:
    #   Tg_eff = (1-eps_f)*Tg + eps_f*Tc_in
    #
    # This is NOT Farokhi-correlations; it’s a class-friendly proxy that preserves monotonicity.

    Kfilm = float(cooling.get('Kfilm', 12.0))  # tunable; higher = more effective for given fcool
    eta_penalty_k = float(cooling.get('eta_penalty_k', 2.0))  # penalty strength

    def wall_temps_for_f(fcool):
        if fcool < 0.0:
            fcool = 0.0
        eps_f = 1.0 - math.exp(-Kfilm * fcool)
        Tg_eff = (1.0 - eps_f) * Tg + eps_f * Tc_in

        # Hot-side convective
        # Assume coolant-side keeps inner wall near Tc_in + small rise; we keep it simple:
        Twc = Tc_in + 0.15 * (Tg_eff - Tc_in)  # inner wall sits closer to coolant than gas

        # conduction
        # Twg = Twc + q"*tw/kw, and q" = hg*(Tg_eff - Twg)
        # => Twg = Twc + hg*(Tg_eff - Twg)*tw/kw
        # => Twg*(1 + hg*tw/kw) = Twc + hg*tw/kw*Tg_eff
        denom = 1.0 + hg * tw / max(1e-12, kw)
        Twg = (Twc + (hg * tw / max(1e-12, kw)) * Tg_eff) / denom

        # Coolant "bulk" exit temperature proxy: rises with extracted heat ~ (Tg_eff - Tc_in)*fcool_scale
        Tc = Tc_in + 0.35 * (Tg_eff - Tc_in)  # conservative bulk
        return Twg, Tc

    # Search minimal fcool satisfying constraints
    fcool = 0.0
    ok = False
    Twg = 1e9
    Tc = 1e9

    # coarse-to-fine
    for step in [0.02, 0.005, 0.001]:
        if fcool < 0.0:
            fcool = 0.0
        while fcool <= 0.30:  # cap coolant fraction
            Twg, Tc = wall_temps_for_f(fcool)
            if (Twg <= Twg_max) and (Tc <= Tc_max):
                ok = True
                break
            fcool += step
        if ok:
            # back up one step and refine
            fcool = max(0.0, fcool - step)
        else:
            break

    # Final evaluate at fcool
    Twg, Tc = wall_temps_for_f(fcool)

    if not ok:
        return fcool, eta_t0 * 0.85, False, "cooling_failed: Twg=%.1f Tc=%.1f" % (Twg, Tc)

    # Apply turbine efficiency penalty for cooling
    # Effective efficiency reduces with coolant fraction: eta_eff = eta0 * (1 - k*fcool)
    eta_t_eff = eta_t0 * max(0.50, 1.0 - eta_penalty_k * fcool)

    return fcool, eta_t_eff, True, "ok: Twg=%.1f Tc=%.1f fcool=%.4f" % (Twg, Tc, fcool)
