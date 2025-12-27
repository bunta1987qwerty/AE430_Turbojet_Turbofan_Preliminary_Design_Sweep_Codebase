# -*- coding: utf-8 -*-
from __future__ import division
import math

from ae430.thermo import R_from_cp_gamma, static_from_total_and_M

def flow_area_from_mdot_total_M(mdot, Tt, pt, M, gamma, cp):
    """
    Returns (A, rho, V, T, p) using isentropic relations.
    Signature matches what your engine.py grep comment expects:
      (mdot, Tt, pt, M, gamma, cp)
    """
    R = R_from_cp_gamma(cp, gamma)
    T, p, rho, a, V = static_from_total_and_M(Tt, pt, M, gamma, R)
    mdot = float(mdot)
    if rho * V <= 1e-12:
        return (float('inf'), rho, V, T, p)
    A = mdot / (rho * V)
    return (A, rho, V, T, p)

def annulus_radii_from_area(A, hub_tip_ratio):
    lam = float(hub_tip_ratio)
    lam = max(1e-6, min(0.999999, lam))
    A = float(A)
    rt = math.sqrt(A / (math.pi * (1.0 - lam * lam)))
    rh = lam * rt
    rm = 0.5 * (rh + rt)
    return rh, rt, rm

def conical_length_from_radii(r1, r2, half_angle_deg):
    theta = math.radians(float(half_angle_deg))
    t = math.tan(theta)
    if abs(t) < 1e-12:
        return float('inf')
    return abs(float(r2) - float(r1)) / t

# Backward-compatible alias some of your patched engine.py tried to import
def conical_diffuser_length(r1, r2, half_angle_deg):
    return conical_length_from_radii(r1, r2, half_angle_deg)
