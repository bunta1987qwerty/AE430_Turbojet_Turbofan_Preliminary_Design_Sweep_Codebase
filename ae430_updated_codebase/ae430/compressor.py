# -*- coding: utf-8 -*-
"""
Compressor module (Python 2.7 compatible).

Provides:
  - compressor_overall(): lumped compressor thermodynamics (Tt, pt, work)
  - estimate_stage_count(): pick stage count based on max per-stage pressure ratio
  - stage_checks_free_vortex(): lightweight geometric/mechanical checks for constraints

NOTE (important):
This project is an engine-cycle + constraints sweep; the compressor "stage checks"
here are deliberately simple / conservative, intended to:
  - flag obviously impossible designs (tip Mach too high, diameter too large, stress too high)
  - provide consistent checks across sweeps
They are not a full mean-line design tool.
"""

from __future__ import division
import math


def compressor_overall(Tt_in, pt_in, pi_c, eta_c, gamma, cp):
    """
    Lumped compressor:
      Tt_out = Tt_in * [ 1 + (pi_c^((gamma-1)/gamma) - 1)/eta_c ]
      pt_out = pt_in * pi_c
      w = cp * (Tt_out - Tt_in)

    Returns dict with Tt_out, pt_out, w, tau_c.
    """
    Tt_in = float(Tt_in)
    pt_in = float(pt_in)
    pi_c = float(pi_c)
    eta_c = float(eta_c)
    gamma = float(gamma)
    cp = float(cp)

    if pi_c <= 1.0:
        # no compression
        return {'Tt_out': Tt_in, 'pt_out': pt_in, 'w': 0.0, 'tau_c': 1.0}

    expn = (gamma - 1.0) / gamma
    tau_is = pi_c ** expn
    tau_c = 1.0 + (tau_is - 1.0) / max(1e-9, eta_c)

    Tt_out = Tt_in * tau_c
    pt_out = pt_in * pi_c
    w = cp * (Tt_out - Tt_in)

    return {'Tt_out': Tt_out, 'pt_out': pt_out, 'w': w, 'tau_c': tau_c}


def estimate_stage_count(pi_c, pi_stage_max):
    """
    Choose smallest integer N such that (pi_c)^(1/N) <= pi_stage_max.
    """
    pi_c = float(pi_c)
    pi_stage_max = float(pi_stage_max)
    if pi_c <= 1.0:
        return 1
    if pi_stage_max <= 1.0:
        return 1
    N = 1
    while True:
        pi_stage = pi_c ** (1.0 / float(N))
        if pi_stage <= pi_stage_max:
            return N
        N += 1
        if N > 200:
            return N


def _speed_of_sound_from_Tt(Tt, gamma, cp):
    """
    Approximate a from Tt using R = cp*(gamma-1)/gamma.
    Using Tt instead of static T is conservative for tip Mach checks.
    """
    Tt = float(Tt)
    gamma = float(gamma)
    cp = float(cp)
    R = cp * (gamma - 1.0) / max(1e-12, gamma)
    return math.sqrt(max(0.0, gamma * R * Tt))


def stage_checks_free_vortex(n_stages,
                            pi_c,
                            Tt_in, pt_in,
                            omega,
                            r_hub, r_tip,
                            Cz,
                            eta_c,
                            gamma, cp,
                            rho_blade,
                            D_max=None,
                            tip_M_max=None,
                            sigma_c_max_MPa=None):
    """
    Lightweight compressor feasibility checks.

    Inputs:
      n_stages: int
      pi_c: overall pressure ratio
      omega: rad/s (your config uses "omega_comp" - treat as rad/s consistently)
      r_hub, r_tip: meters
      Cz: axial velocity proxy (m/s) - included for future extension; not used heavily here
      rho_blade: kg/m^3
      D_max: max allowed tip diameter (m)
      tip_M_max: max allowed blade tip Mach number
      sigma_c_max_MPa: max allowed blade stress proxy (MPa)

    Returns dict:
      ok: bool
      max_Mtip: float
      sigma_MPa: float
      D_tip: float
      pi_stage: float
      notes: list[str]
    """
    notes = []
    n_stages = int(n_stages)
    pi_c = float(pi_c)
    Tt_in = float(Tt_in)
    omega = float(omega)
    r_hub = float(r_hub)
    r_tip = float(r_tip)
    rho_blade = float(rho_blade)

    if r_tip <= 0.0 or r_tip <= r_hub:
        return {
            'ok': False,
            'max_Mtip': float('nan'),
            'sigma_MPa': float('nan'),
            'D_tip': float('nan'),
            'pi_stage': float('nan'),
            'notes': ['bad_radii']
        }

    # Stage pressure ratio check
    if n_stages <= 0:
        n_stages = 1
    pi_stage = 1.0
    if pi_c > 1.0:
        pi_stage = pi_c ** (1.0 / float(n_stages))

    # Diameter check
    D_tip = 2.0 * r_tip
    if D_max is not None:
        if D_tip > float(D_max) + 1e-12:
            notes.append('D_tip_exceeds_D_max')

    # Tip Mach check
    U_tip = omega * r_tip  # m/s
    a_tip = _speed_of_sound_from_Tt(Tt_in, gamma, cp)
    M_tip = U_tip / max(1e-12, a_tip)

    if tip_M_max is not None:
        if M_tip > float(tip_M_max) + 1e-12:
            notes.append('tip_M_exceeds_limit')

    # Stress proxy (very simple):
    # Use dynamic-pressure-like scaling: sigma ~ 0.5*rho*U^2
    # Convert Pa -> MPa
    sigma_Pa = 0.5 * rho_blade * (U_tip ** 2)
    sigma_MPa = sigma_Pa / 1e6

    if sigma_c_max_MPa is not None:
        if sigma_MPa > float(sigma_c_max_MPa) + 1e-12:
            notes.append('sigma_exceeds_limit')

    ok = (len(notes) == 0)

    return {
        'ok': ok,
        'max_Mtip': M_tip,
        'sigma_MPa': sigma_MPa,
        'D_tip': D_tip,
        'pi_stage': pi_stage,
        'notes': notes
    }
