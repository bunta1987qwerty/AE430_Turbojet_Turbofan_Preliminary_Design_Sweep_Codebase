# -*- coding: utf-8 -*-
"""
Inlet / diffuser + installed terms (Python 2.7).

Implements project-installed thrust relation:
  F_inst = F_uninst - D_add + F_lip

We return:
  pi_d (diffuser pressure recovery),
  D_add, F_lip,
  and notes.

To keep this compatible with the rest of the pipeline, we expose:
  inlet_performance(M0, q0, A_capture, cfg_inlet, AR=None, phi_deg=None)

You can drive diffuser recovery by sweeping inlet.pi_d_factor (constant model),
or by setting inlet.model="sk_proxy" and passing AR, phi_deg from geometry.
"""

from __future__ import division
from ae430.diffuser import diffuser_pressure_recovery

def inlet_performance(M0, q0, A_capture, cfg_inlet, AR=None, phi_deg=None):
    notes = []
    if cfg_inlet is None:
        cfg_inlet = {}

    Cd_add = float(cfg_inlet.get('Cd_add', 0.0))
    Cl_lip = float(cfg_inlet.get('Cl_lip', 0.0))

    D_add = Cd_add * float(q0) * float(A_capture)
    F_lip = Cl_lip * float(q0) * float(A_capture)

    pi_d = diffuser_pressure_recovery(cfg_inlet, AR=AR, phi_deg=phi_deg,
                                      fallback_pi=float(cfg_inlet.get('pi_d', 0.98)))

    notes.append("inlet.model=%s pi_d=%.6f" % (str(cfg_inlet.get('model', 'constant')), pi_d))
    return {'pi_d': pi_d, 'D_add': D_add, 'F_lip': F_lip, 'notes': notes}
