# -*- coding: utf-8 -*-
from __future__ import division

import copy
import math
import os

from ae430.thermo import speed_of_sound, total_T_from_static, total_p_from_static, static_from_total_and_M
from ae430.geometry import flow_area_from_mdot_total_M, annulus_radii_from_area
from ae430.inlet import inlet_performance
from ae430.diffuser import diffuser_pressure_recovery
from ae430.compressor import compressor_overall, estimate_stage_count, stage_checks_free_vortex
from ae430.cantera_burner import find_max_f_by_constraints, equilibrium_burner_state
from ae430.cooling import solve_cooling_and_eta
from ae430.turbine import turbine_from_required_power
from ae430.nozzle import nozzle_expand_to_ambient
from ae430.cost import cost_function


def _nan():
    try:
        return float('nan')
    except Exception:
        return 0.0 / 0.0


def _push_reason(reasons, s):
    if s is None:
        return
    reasons.append(str(s))


def _resolve_path(cfg, maybe_rel_path):
    """
    Resolve files like diesel.xml relative to cfg['_base_dir'] if present.
    """
    if maybe_rel_path is None:
        return None
    p = str(maybe_rel_path)
    if os.path.isabs(p):
        return p
    base_dir = cfg.get('_base_dir', None)
    if base_dir:
        cand = os.path.join(base_dir, p)
        if os.path.exists(cand):
            return cand
    # fallback: cwd
    return p


def _compute_freestream(cfg):
    cst = cfg.get('constants', {})
    gamma_cold = float(cst.get('gamma_cold', 1.4))
    cp_cold = float(cst.get('cp_cold', 1000.0))
    R_cold = cp_cold * (gamma_cold - 1.0) / gamma_cold

    fl = cfg.get('flight', {})
    T0 = float(fl.get('T0', 220.0))
    p0 = float(fl.get('p0', 23000.0))
    V0 = float(fl.get('V0', 250.0))

    a0 = speed_of_sound(T0, gamma_cold, R_cold)
    M0 = V0 / a0 if a0 > 1e-12 else 0.0

    Tt0 = total_T_from_static(T0, M0, gamma_cold)
    pt0 = total_p_from_static(p0, M0, gamma_cold)

    # dynamic pressure proxy q = 0.5*rho V^2, using q = 0.5*gamma*p*M^2 (ideal gas)
    q0 = 0.5 * gamma_cold * p0 * (M0 ** 2)

    return {
        'gamma_cold': gamma_cold, 'cp_cold': cp_cold, 'R_cold': R_cold,
        'T0': T0, 'p0': p0, 'V0': V0, 'a0': a0, 'M0': M0, 'Tt0': Tt0, 'pt0': pt0, 'q0': q0
    }


def _cycle_given_mdot(cfg, mdot0, logger):
    """
    Run one full cycle for a given mdot0.
    Returns: (results_dict, stations_dict, constraints_dict)
    """
    reasons = []

    req = cfg.get('requirements', {})
    losses = cfg.get('losses', {})
    design = cfg.get('design', {})
    geom = cfg.get('geometry', {})
    asm = cfg.get('assumptions', {})
    cst = cfg.get('constants', {})

    fs = _compute_freestream(cfg)
    gamma_cold = fs['gamma_cold']
    cp_cold = fs['cp_cold']
    R_cold = fs['R_cold']

    gamma_hot = float(cst.get('gamma_hot', 1.3))
    cp_hot = float(cst.get('cp_hot', 1250.0))
    R_hot = cp_hot * (gamma_hot - 1.0) / gamma_hot

    engine_type = design.get('engine_type', 'turbojet')
    alpha = float(design.get('alpha', 0.0))
    if engine_type != 'turbofan':
        alpha = 0.0

    # sweepable
    pi_c = float(design.get('pi_c', 25.0))
    pi_f = float(design.get('pi_f', 1.6))

    eta_c = float(losses.get('eta_c', 0.90))
    eta_f = float(losses.get('eta_f', 0.90))
    eta_m = float(losses.get('eta_m', 0.98))
    eta_n = float(losses.get('eta_n', 0.95))
    pi_b = float(losses.get('pi_b', 0.98))

    # spool / structural constraint inputs
    omega_comp = float(cfg.get('spools', {}).get('omega_comp', 500.0))
    rho_blade = float(asm.get('rho_blade', 8000.0))
    D_max = float(cfg.get('compressor', {}).get('D_max', 0.60))
    pi_stage_max = float(cfg.get('compressor', {}).get('pi_stage_max', 1.25))
    tip_M_max = float(req.get('tip_M_max', 1.2))
    sigma_c_max = float(req.get('sigma_c_max_MPa', 70.0))

    # inlet / capture
    inlet_cfg = cfg.get('inlet', {})
    M2 = float(asm.get('M2', 0.4))
    phi_in_deg = float(geom.get('inlet_diffuser_half_angle_deg', 7.0))

    # Capture area based on station-0 totals and assumed M2 at face
    A_capture, _, _ = flow_area_from_mdot_total_M(mdot0, fs['Tt0'], fs['pt0'], M2, gamma_cold, cp_cold)

    # iterate to compute pi_d based on AR
    pi_d = float(inlet_cfg.get('pi_d', 0.98))
    Tt2 = fs['Tt0']
    pt2 = fs['pt0'] * pi_d

    inlet_out = None
    for _ in range(3):
        A2, _, _ = flow_area_from_mdot_total_M(mdot0, Tt2, pt2, M2, gamma_cold, cp_cold)
        AR_in = A2 / max(1e-12, A_capture)
        inlet_out = inlet_performance(fs['M0'], fs['q0'], A_capture, inlet_cfg, AR=AR_in, phi_deg=phi_in_deg)
        pi_d = inlet_out.get('pi_d', pi_d)
        pt2 = fs['pt0'] * pi_d

    if inlet_out is None:
        inlet_out = {'pi_d': pi_d, 'D_add': 0.0, 'F_lip': 0.0}

    D_add = float(inlet_out.get('D_add', 0.0))
    F_lip = float(inlet_out.get('F_lip', 0.0))

    st = {}
    st['0'] = {'Tt': fs['Tt0'], 'pt': fs['pt0']}
    st['2'] = {'Tt': Tt2, 'pt': pt2}

    # split
    mdot_core = mdot0 / (1.0 + alpha)
    mdot_byp = mdot0 - mdot_core

    # fan for turbofan
    fan = None
    if alpha > 0.0:
        fan = compressor_overall(st['2']['Tt'], st['2']['pt'], pi_f, eta_f, gamma_cold, cp_cold)
        st['13'] = {'Tt': fan['Tt_out'], 'pt': fan['pt_out']}
    else:
        st['13'] = {'Tt': st['2']['Tt'], 'pt': st['2']['pt']}

    # core compressor
    comp = compressor_overall(st['2']['Tt'], st['2']['pt'], pi_c, eta_c, gamma_cold, cp_cold)
    st['3'] = {'Tt': comp['Tt_out'], 'pt': comp['pt_out']}

    # stage count and checks
    n_stages = cfg.get('compressor', {}).get('n_stages', None)
    if n_stages is None:
        n_stages = estimate_stage_count(pi_c, pi_stage_max)
    n_stages = int(n_stages)

    # geometry radii for checks (use A2 and hub-tip)
    lam_comp = float(geom.get('hub_tip_ratio_comp', 0.6))
    A2, _, _ = flow_area_from_mdot_total_M(mdot0, st['2']['Tt'], st['2']['pt'], M2, gamma_cold, cp_cold)
    rh2, rt2, _ = annulus_radii_from_area(A2, lam_comp)

    Cz_comp = float(asm.get('Cz_comp', 150.0))
    comp_checks = stage_checks_free_vortex(
        n_stages=n_stages, pi_c=pi_c,
        Tt_in=st['2']['Tt'], pt_in=st['2']['pt'],
        omega=omega_comp,
        r_hub=rh2, r_tip=rt2, Cz=Cz_comp,
        eta_c=eta_c, gamma=gamma_cold, cp=cp_cold,
        rho_blade=rho_blade,
        D_max=D_max,
        tip_M_max=tip_M_max,
        sigma_c_max_MPa=sigma_c_max
    )

    # burner diffuser (3 -> 3b)
    bd_cfg = cfg.get('burner_diffuser', {})
    phi_bd_deg = float(geom.get('burner_diffuser_half_angle_deg', 7.0))
    M3 = float(asm.get('M3', 0.4))
    M_b = float(asm.get('M_burner', 0.1))

    A3, _, _ = flow_area_from_mdot_total_M(mdot_core, st['3']['Tt'], st['3']['pt'], M3, gamma_cold, cp_cold)
    # first guess
    pi_bd = diffuser_pressure_recovery(bd_cfg, AR=1.0, phi_deg=phi_bd_deg, fallback_pi=float(bd_cfg.get('pi_bd', 0.99)))
    pt_bd = st['3']['pt'] * pi_bd
    Ab, _, _ = flow_area_from_mdot_total_M(mdot_core, st['3']['Tt'], pt_bd, M_b, gamma_cold, cp_cold)
    AR_bd = Ab / max(1e-12, A3)
    pi_bd = diffuser_pressure_recovery(bd_cfg, AR=AR_bd, phi_deg=phi_bd_deg, fallback_pi=float(bd_cfg.get('pi_bd', 0.99)))
    pt_bd = st['3']['pt'] * pi_bd
    st['3b'] = {'Tt': st['3']['Tt'], 'pt': pt_bd}

    # burner via cantera equilibrium search
    burn = cfg.get('burner', {})
    f_search = burn.get('f_search', {'f_min': 0.0, 'f_max': 0.08, 'n_grid': 60})
    Tt4_target = float(burn.get('Tt4_target', 1850.0))
    Tt4_max = float(req.get('Tt4_max', burn.get('Tt4_max', 1900.0)))
    XCO_max = float(req.get('XCO_max', burn.get('XCO_max', 2e-4)))

    mech_xml = _resolve_path(cfg, burn.get('mech_xml', 'diesel.xml'))
    fuel_species = burn.get('fuel_species', 'C12H24-1')

    f_res = find_max_f_by_constraints(
        Tt3=st['3b']['Tt'],
        pt3=st['3b']['pt'],
        pi_b=pi_b,
        mech_xml_path=mech_xml,
        fuel_species=fuel_species,
        f_min=float(f_search.get('f_min', 0.0)),
        f_max=float(f_search.get('f_max', 0.08)),
        n_grid=int(f_search.get('n_grid', 60)),
        Tt4_target=Tt4_target,
        Tt4_max=Tt4_max,
        XCO_max=XCO_max
    )

    if (f_res is None) or (not bool(f_res.get('ok', False))):
        _push_reason(reasons, "burner_no_feasible_f")
        results = {'feasible': False, 'fail_reasons': reasons, 'mdot0': mdot0}
        cons = _constraints_bundle(cfg, inlet_out, comp_checks, None, None, None, None)
        return results, st, cons

    f = float(f_res.get('f', 0.0))
    burn_state = equilibrium_burner_state(
        Tt3=st['3b']['Tt'], pt3=st['3b']['pt'],
        f=f, pi_b=pi_b,
        mech_xml_path=mech_xml,
        fuel_species=fuel_species
    )
    st['4'] = {'Tt': burn_state['Tt4'], 'pt': burn_state['pt4']}

    # cooling
    fcool, eta_t_eff, cool_ok, cool_reason = solve_cooling_and_eta(cfg, st['4']['Tt'], st['3']['Tt'])

    # shaft work requirement (per kg core)
    w_comp = float(comp.get('w', 0.0))  # J/kg
    w_fan = 0.0
    if alpha > 0.0 and fan is not None:
        w_fan = float(fan.get('w', 0.0)) * (mdot_byp / max(1e-12, mdot_core))
    w_req = (w_comp + w_fan) / max(1e-12, eta_m)

    # turbine match
    turb = turbine_from_required_power(
        Tt4=st['4']['Tt'], pt4=st['4']['pt'],
        w_required=w_req,
        f=f,
        eta_t=eta_t_eff,
        gamma_hot=gamma_hot,
        cp_hot=cp_hot
    )
    if not bool(turb.get('ok', False)):
        _push_reason(reasons, turb.get('fail_reason', 'turbine_fail'))
        results = {'feasible': False, 'fail_reasons': reasons, 'mdot0': mdot0}
        cons = _constraints_bundle(cfg, inlet_out, comp_checks, f_res, burn_state, (fcool, cool_ok, cool_reason), turb)
        return results, st, cons

    st['5'] = {'Tt': turb['Tt5'], 'pt': turb['pt5']}

    # nozzles
    p_amb = fs['p0']
    nozzle_half_angle = float(geom.get('nozzle_half_angle_deg', 7.0))

    noz_core = nozzle_expand_to_ambient(
        pt_in=st['5']['pt'], Tt_in=st['5']['Tt'], p_amb=p_amb,
        gamma=gamma_hot, cp=cp_hot, eta_n=eta_n,
        mdot=mdot_core * (1.0 + f),
        half_angle_deg=nozzle_half_angle,
        assume_perfect_expansion=True
    )
    noz_byp = None
    if alpha > 0.0:
        noz_byp = nozzle_expand_to_ambient(
            pt_in=st['13']['pt'], Tt_in=st['13']['Tt'], p_amb=p_amb,
            gamma=gamma_cold, cp=cp_cold, eta_n=eta_n,
            mdot=mdot_byp,
            half_angle_deg=nozzle_half_angle,
            assume_perfect_expansion=True
        )

    V0 = fs['V0']
    F_core = float(noz_core.get('F', 0.0))
    F_byp = 0.0 if noz_byp is None else float(noz_byp.get('F', 0.0))

    # uninstalled: (gross thrust) - mdot0*V0
    F_uninst = (F_core + F_byp) - mdot0 * V0
    F_inst = F_uninst - D_add + F_lip

    # residence time constraint (very simple burner length proxy)
    tau_req = float(req.get('tau_residence_min', 0.005))
    Ts_b, ps_b, rho_b, V_b = static_from_total_and_M(st['3b']['Tt'], st['3b']['pt'], M_b, gamma_cold, cp_cold)
    lb_min = float(geom.get('lb_min', 0.2))
    lam_b = float(geom.get('hub_tip_ratio_burner', 0.7))
    rhb, rtb, _ = annulus_radii_from_area(Ab, lam_b)
    D_b = 2.0 * rtb
    spacing = float(geom.get('spacing_factor', 0.25))
    L_b = max(lb_min, spacing * D_b)
    tau_res = L_b / max(1e-12, V_b)

    # overall efficiency proxy (propulsive power / fuel power)
    LHV = float(cfg.get('fuel', {}).get('LHV', 43e6))
    mdot_f = mdot_core * f
    eta0 = 0.0
    if mdot_f > 1e-12:
        eta0 = (F_inst * V0) / (mdot_f * LHV)

    # cost
    Amax = max(A_capture, A2, A3, Ab)
    Nspools = int(design.get('Nspools', 2))
    fc = cost_function(Amax, 1.0, eta0, Nspools)  # length proxy fixed; your cost.py already handles

    cons = _constraints_bundle(cfg, inlet_out, comp_checks, f_res, burn_state, (fcool, cool_ok, cool_reason), turb)
    cons['tau_residence_min'] = {'ok': (tau_res >= tau_req), 'value': tau_res, 'limit': tau_req}

    feasible = True
    for k in cons.keys():
        if not bool(cons[k].get('ok', True)):
            feasible = False

    if not feasible:
        _push_reason(reasons, "constraint_violation")

    results = {
        'feasible': feasible,
        'fail_reasons': reasons,
        'engine_type': engine_type,
        'alpha': alpha,
        'mdot0': mdot0,
        'mdot_core': mdot_core,
        'mdot_bypass': mdot_byp,
        'pi_d': pi_d,
        'pi_bd': pi_bd,
        'pi_c': pi_c,
        'pi_f': pi_f,
        'f': f,
        'fcool': fcool,
        'eta_t_eff': eta_t_eff,
        'pi_t': turb.get('pi_t', _nan()),
        'Tt4': burn_state.get('Tt4', _nan()),
        'XCO': burn_state.get('X_CO', burn_state.get('XCO', _nan())),
        'F_uninstalled': F_uninst,
        'F_installed': F_inst,
        'eta0': eta0,
        'tau_res': tau_res,
        'A_capture': A_capture,
        'A2': A2,
        'A3': A3,
        'Ab': Ab,
        'Amax_engine': Amax,
        'fc': fc,
    }
    return results, st, cons


def _constraints_bundle(cfg, inlet_out, comp_checks, f_res, burn_state, cooling_tuple, turb):
    req = cfg.get('requirements', {})
    cons = {}

    # compressor checks
    cons['tip_M_max'] = {
        'ok': bool(comp_checks.get('max_Mtip', 0.0) <= float(req.get('tip_M_max', 1.2))),
        'value': comp_checks.get('max_Mtip', _nan()),
        'limit': float(req.get('tip_M_max', 1.2))
    }
    cons['sigma_c_max_MPa'] = {
        'ok': bool(comp_checks.get('sigma_MPa', 0.0) <= float(req.get('sigma_c_max_MPa', 70.0))),
        'value': comp_checks.get('sigma_MPa', _nan()),
        'limit': float(req.get('sigma_c_max_MPa', 70.0))
    }
    cons['compressor_stage_checks'] = {
        'ok': bool(comp_checks.get('ok', False)),
        'value': 1 if comp_checks.get('ok', False) else 0,
        'limit': 1,
        'notes': "; ".join(comp_checks.get('notes', []))
    }

    # burner
    if burn_state is not None:
        cons['Tt4_max'] = {
            'ok': bool(float(burn_state.get('Tt4', 1e9)) <= float(req.get('Tt4_max', 1900.0))),
            'value': float(burn_state.get('Tt4', _nan())),
            'limit': float(req.get('Tt4_max', 1900.0))
        }
        XCO = burn_state.get('X_CO', burn_state.get('XCO', _nan()))
        cons['XCO_max'] = {
            'ok': bool(float(XCO) <= float(req.get('XCO_max', 2e-4))),
            'value': float(XCO),
            'limit': float(req.get('XCO_max', 2e-4))
        }
    else:
        cons['Tt4_max'] = {'ok': False, 'value': _nan(), 'limit': float(req.get('Tt4_max', 1900.0))}
        cons['XCO_max'] = {'ok': False, 'value': _nan(), 'limit': float(req.get('XCO_max', 2e-4))}

    # cooling
    if cooling_tuple is not None:
        fcool, cool_ok, cool_reason = cooling_tuple
        cons['cooling_ok'] = {'ok': bool(cool_ok), 'value': 1 if cool_ok else 0, 'limit': 1, 'notes': cool_reason}
    else:
        cons['cooling_ok'] = {'ok': False, 'value': 0, 'limit': 1, 'notes': 'missing'}

    # turbine ok
    if turb is not None:
        cons['turbine_ok'] = {
            'ok': bool(turb.get('ok', False)),
            'value': 1 if turb.get('ok', False) else 0,
            'limit': 1,
            'notes': turb.get('fail_reason', '')
        }
    else:
        cons['turbine_ok'] = {'ok': False, 'value': 0, 'limit': 1, 'notes': 'missing'}

    return cons


def evaluate_design(cfg, logger):
    """
    Solve mdot to match required installed thrust.

    Returns (results, stations, constraints).
    """
    F_required = float(cfg.get('requirements', {}).get('F_required', 350e3))

    # robust secant with fallback
    md0 = float(cfg.get('design', {}).get('mdot_guess', 325.0))
    md1 = md0 * 1.10

    def eval_F(mdot):
        res, st, cons = _cycle_given_mdot(cfg, mdot, logger)
        if not bool(res.get('feasible', False)) and ('burner_no_feasible_f' in res.get('fail_reasons', [])):
            return -1e12, res, st, cons
        return float(res.get('F_installed', -1e12)), res, st, cons

    F0, r0, s0, c0 = eval_F(md0)
    F1, r1, s1, c1 = eval_F(md1)

    # if both completely dead
    if (F0 < -1e11) and (F1 < -1e11):
        r0['feasible'] = False
        _push_reason(r0.setdefault('fail_reasons', []), "F_solve_failed")
        return r0, s0, c0

    best = (abs(F0 - F_required), r0, s0, c0)
    if abs(F1 - F_required) < best[0]:
        best = (abs(F1 - F_required), r1, s1, c1)

    for _ in range(12):
        if abs(F1 - F0) < 1e-9:
            break

        md2 = md1 + (F_required - F1) * (md1 - md0) / (F1 - F0)
        if md2 < 1.0:
            md2 = 1.0
        if md2 > 3000.0:
            md2 = 3000.0

        F2, r2, s2, c2 = eval_F(md2)
        if abs(F2 - F_required) < best[0]:
            best = (abs(F2 - F_required), r2, s2, c2)

        # convergence
        if abs(F2 - F_required) / max(1.0, F_required) < 1e-3:
            r2['F_required'] = F_required
            return r2, s2, c2

        md0, F0 = md1, F1
        md1, F1 = md2, F2

    # return best attempt but mark if thrust mismatch big
    _, rb, sb, cb = best
    rb['F_required'] = F_required
    if abs(float(rb.get('F_installed', -1e12)) - F_required) / max(1.0, F_required) > 0.05:
        rb['feasible'] = False
        _push_reason(rb.setdefault('fail_reasons', []), "F_mismatch")
    return rb, sb, cb
