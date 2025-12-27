# -*- coding: utf-8 -*-
"""
Diffuser pressure recovery helper.

We need a way to model total-pressure recovery (pi = pt_out/pt_in) for:
  - inlet diffuser  (station 0 -> 2)
  - burner diffuser (station 3 -> 3b)

Class materials often show a chart/correlation (digitized). If you do not
have a digitized correlation, the sweep approach is acceptable: sweep a
factor pi_*_factor within a reasonable range, and clamp to min/max.

This module supports:
  - model='constant': fixed pi from JSON (with optional factor + clamp)
  - model='sk_proxy': smooth proxy loss model that depends on AR and half-angle

All outputs are clamped to [pi_min, pi_max].
Python 2.7 compatible.
"""
from __future__ import division
import math

def _clamp(x, lo, hi):
    if x < lo:
        return lo
    if x > hi:
        return hi
    return x

def diffuser_pressure_recovery(cfg_section, AR=None, half_angle_rad=None, phi_deg=None, fallback_pi=0.98):
    """
    Parameters
    ----------
    cfg_section : dict
        e.g. cfg['inlet'] or cfg['burner_diffuser']
        expected fields (depending on model):
          model: 'constant' or 'sk_proxy'
          pi_d or pi_bd (base)
          pi_d_factor or pi_bd_factor (multiplier)
          pi_*_min, pi_*_max
          loss_min, k_ar, k_phi, k_short, phi_crit_deg

    AR : float or None
        Area ratio A_out/A_in (>=1 typical for diffusers)
    half_angle_rad : float or None
        Diffuser half-angle in radians (optional)
    phi_deg : float or None
        Diffuser half-angle in degrees (optional)

    Returns
    -------
    pi : float
        Pressure recovery pt_out/pt_in
    """
    if cfg_section is None:
        cfg_section = {}

    model = cfg_section.get('model', 'constant')

    # Choose base + factor keys depending on whether this is inlet or burner diffuser
    # (both configs exist in your JSON as 'pi_d' and 'pi_bd')
    base = None
    factor = 1.0
    pi_min = None
    pi_max = None

    if 'pi_d' in cfg_section or 'pi_d_min' in cfg_section or 'pi_d_max' in cfg_section:
        base = float(cfg_section.get('pi_d', fallback_pi))
        factor = float(cfg_section.get('pi_d_factor', cfg_section.get('pi_d_factor', 1.0)))
        pi_min = float(cfg_section.get('pi_d_min', 0.0))
        pi_max = float(cfg_section.get('pi_d_max', 1.0))
    elif 'pi_bd' in cfg_section or 'pi_bd_min' in cfg_section or 'pi_bd_max' in cfg_section:
        base = float(cfg_section.get('pi_bd', fallback_pi))
        factor = float(cfg_section.get('pi_bd_factor', cfg_section.get('pi_bd_factor', 1.0)))
        pi_min = float(cfg_section.get('pi_bd_min', 0.0))
        pi_max = float(cfg_section.get('pi_bd_max', 1.0))
    else:
        base = float(cfg_section.get('pi', fallback_pi))
        factor = float(cfg_section.get('pi_factor', 1.0))
        pi_min = float(cfg_section.get('pi_min', 0.0))
        pi_max = float(cfg_section.get('pi_max', 1.0))

    # Convert angle to degrees if needed
    if phi_deg is None and half_angle_rad is not None:
        phi_deg = float(half_angle_rad) * 180.0 / math.pi
    if phi_deg is None:
        phi_deg = float(cfg_section.get('phi_deg', cfg_section.get('phi_crit_deg', 0.0)))

    if model == 'constant':
        pi = base * factor
        return _clamp(pi, pi_min, pi_max)

    if model == 'sk_proxy':
        # Smooth proxy loss model:
        # loss = loss_min + k_ar*(AR-1)^2 + k_phi*max(0,phi-phi_crit)^2 + k_short/AR
        loss_min = float(cfg_section.get('loss_min', 0.01))
        k_ar = float(cfg_section.get('k_ar', 0.02))
        k_phi = float(cfg_section.get('k_phi', 0.02))
        k_short = float(cfg_section.get('k_short', 0.01))
        phi_crit = float(cfg_section.get('phi_crit_deg', 6.0))

        if AR is None:
            AR = 1.0
        AR = float(AR)
        if AR < 1.0:
            # treat contraction as harsher
            AR_eff = 1.0 / max(1e-9, AR)
        else:
            AR_eff = AR

        dAR = (AR_eff - 1.0)
        dphi = phi_deg - phi_crit
        if dphi < 0.0:
            dphi = 0.0

        loss = loss_min + k_ar * (dAR ** 2) + k_phi * (dphi ** 2) + (k_short / max(1.0, AR_eff))
        pi = (1.0 - loss) * factor

        return _clamp(pi, pi_min, pi_max)

    # Unknown model
    pi = base * factor
    return _clamp(pi, pi_min, pi_max)
