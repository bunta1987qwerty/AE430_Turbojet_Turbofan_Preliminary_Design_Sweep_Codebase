# -*- coding: utf-8 -*-
from __future__ import division
import math

def R_from_cp_gamma(cp, gamma):
    return cp * (gamma - 1.0) / gamma

def speed_of_sound(T, gamma, R):
    return math.sqrt(max(0.0, gamma * R * T))

def total_T_from_static(T, M, gamma):
    return T * (1.0 + 0.5 * (gamma - 1.0) * M * M)

def total_p_from_static(p, M, gamma):
    return p * (1.0 + 0.5 * (gamma - 1.0) * M * M) ** (gamma / (gamma - 1.0))

def static_T_from_total(Tt, M, gamma):
    return Tt / (1.0 + 0.5 * (gamma - 1.0) * M * M)

def static_p_from_total(pt, Tt, T, gamma):
    # pt/p = (Tt/T)^(gamma/(gamma-1))
    return pt / ((Tt / T) ** (gamma / (gamma - 1.0)))

def static_from_total_and_M(Tt, pt, M, gamma, R):
    T = static_T_from_total(Tt, M, gamma)
    p = static_p_from_total(pt, Tt, T, gamma)
    rho = p / (R * T)
    a = speed_of_sound(T, gamma, R)
    V = M * a
    return T, p, rho, a, V
