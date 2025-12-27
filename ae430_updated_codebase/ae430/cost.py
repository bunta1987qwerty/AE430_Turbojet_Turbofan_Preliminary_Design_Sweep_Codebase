# -*- coding: utf-8 -*-
"""
Cost function and bookkeeping.

Project statement defines:
  - fc as a function of max cross-sectional area, engine length, overall efficiency, and Nspools
  - l = l_id + l_c + l_bd + l_b + l_t + l_n

The PDF formatting can be hard to parse when extracted as text, but the intended form is:

    f_c = sqrt(A_max,engine) * sqrt(l) / (eta_0^2) * (3 + N_spools)

You should verify this matches your handout.
"""
from __future__ import division
import math

def cost_function(Amax_engine, l_total, eta0, Nspools):
    """
    Compute f_c.

    Returns nan if eta0 <= 0.
    """
    if eta0 <= 0.0:
        return float('nan')
    return (math.sqrt(Amax_engine) * math.sqrt(l_total) / (eta0**2)) * (3.0 + float(Nspools))
