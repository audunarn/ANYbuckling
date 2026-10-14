'''
Regression tests for the semi-analytical S3/U3 fixes of the calculation review
(docs/review/CALCULATION_REVIEW.md, findings U1-U9 and S1-S7).

Every test compares with an independent reference -- closed-form plate
buckling coefficients, an independent Ritz solution, the Euler column of the
stiffener with its full attached plate, sigma/fy for tension yield, the
DNV-RP-C201 (2010) and IACS/CSR formulas evaluated by hand, or monotonicity and
invariance properties -- never with the solver's own output.  Each test fails
on the reviewed revision of ``solver.py``.
'''
import dataclasses
import math
from types import SimpleNamespace

import numpy as np
import pytest

from anybuckling import FlatStru
from anybuckling.semianalytical import solver as S

E, NU, FY = 210000.0, 0.3, 355.0
README_BULB = (260, 12, 49, 28, 'bulb')


def flat(s, t, l, stf, sx=0.0, sx2=None, sy=0.0, tau=0.0, p=0.0, method='SemiAnalytical S3/U3',
         stiffener_support='Continuous', **buckling_parameters):
    '''FlatStru stiffened panel, as a user would call it (mm, MPa).'''
    panel = FlatStru('Flat plate, stiffened')
    panel.set_material(FY, E, 1.15, NU)
    panel.set_plate_geometry(spacing=s, thickness=t, span=l)
    panel.set_stresses(pressure=p, sigma_x1=sx, sigma_x2=sx if sx2 is None else sx2, sigma_y1=sy, sigma_y2=sy,
                       tau_xy=tau)
    hw, tw, bf, tf, stf_type = stf
    panel.set_stiffener(hw=hw, tw=tw, bf=bf, tf=tf, stf_type=stf_type, spacing=s)
    panel.set_buckling_parameters(calculation_method=method, buckling_acceptance='ultimate',
                                  stiffener_support=stiffener_support, **buckling_parameters)
    if method.startswith('Semi'):
        panel.set_puls_parameters(stiffener_end=stiffener_support)
    return panel


def u3_row(a, b, t, sx1, sx2=None, sy1=0.0, sy2=None, tau=0.0, p=0.0, rot=('SS', 'SS'), support='Integrated'):
    return {'Plate length': a, 'Plate width': b, 'Plate thick.': t, 'Yield stress plate': FY,
            'Axial stress': sx1, 'Axial stress 2': sx1 if sx2 is None else sx2,
            'Trans. stress 1': sy1, 'Trans. stress 2': sy1 if sy2 is None else sy2, 'Shear stress': tau,
            'Pressure (fixed)': p, 'In-plane support': support,
            'Rotational support': rot[0], 'Rotational support 2': rot[1]}


def s3_panel(**overrides):
    values = dict(length=3300.0, stiffener_spacing=680.0, plate_thickness=12.0, stiffener_type='L-bulb',
                  stiffener_boundary='Cont', stiffener_height=260.0, web_thickness=12.0, flange_width=49.0,
                  flange_thickness=28.0, yield_stress_plate=FY, yield_stress_stiffener=FY, axial_stress=50.0,
                  transverse_stress_1=100.0, transverse_stress_2=100.0, shear_stress=5.0, pressure=0.01,
                  in_plane_support='Integrated')
    values.update(overrides)
    return S.S3PanelInput(**values)


def sigma_e(b, t):
    '''Classical plate reference stress pi^2 D / (b^2 t).'''
    return math.pi ** 2 * E / (12.0 * (1.0 - NU ** 2)) * (t / b) ** 2


def k_uniform(a, b):
    '''Classical SS plate under uniaxial compression: k = min_m (m b/a + a/(m b))^2.'''
    return min((m * b / a + a / (m * b)) ** 2 for m in range(1, 200))


def euler_full_plate(s, t, length, hw, tw, bf, tf):
    '''Euler column stress of the stiffener with its full attached plate (pinned over the span).'''
    parts = [(s * t, 0.0, s * t ** 3 / 12.0), (hw * tw, t / 2 + hw / 2, tw * hw ** 3 / 12.0),
             (bf * tf, t / 2 + hw + tf / 2, bf * tf ** 3 / 12.0)]
    area = sum(p[0] for p in parts)
    centroid = sum(p[0] * p[1] for p in parts) / area
    inertia = sum(i + ar * (z - centroid) ** 2 for ar, z, i in parts)
    return math.pi ** 2 * E * inertia / (area * length ** 2)


def k_gradient_ritz(alpha, psi, mmax=25, nmax=12):
    '''Independent Ritz solution: SS plate, sigma_x varying linearly across the width.

    N_x(y) = N_1 (1 - (1 - psi) y / b); returns k with sigma_1,cr = k pi^2 D / (b^2 t).
    '''
    a, b = alpha, 1.0
    best = None
    for m in range(1, mmax + 1):
        kx = m * math.pi / a
        lin = np.zeros((nmax, nmax))
        geo = np.zeros((nmax, nmax))
        for i in range(nmax):
            n = i + 1
            lin[i, i] = a * b / 4.0 * (kx * kx + (n * math.pi / b) ** 2) ** 2
            for j in range(nmax):
                n2 = j + 1
                if n == n2:
                    integral = b / 2.0 - (1.0 - psi) * b / 4.0
                elif (n + n2) % 2:
                    integral = (1.0 - psi) * b / math.pi ** 2 * 4.0 * n * n2 / (n * n - n2 * n2) ** 2
                else:
                    integral = 0.0
                geo[i, j] = kx * kx * a / 2.0 * integral
        scale = np.diag(1.0 / np.sqrt(np.diag(lin)))
        top = np.linalg.eigvalsh(scale @ geo @ scale).max()
        if top > 0.0:
            k = 1.0 / top / math.pi ** 2
            best = k if best is None else min(best, k)
    return best


def u3_elastic_factor(panel, config=None):
    config = config or S.S3SolverConfig()
    return S.elastic_u3_buckling_factors(panel, S.build_u3_ritz_modes(panel, config), config)['critical_factor']


def _u3_structure(up_boundary='SSSS', **extra):
    plate = SimpleNamespace(
        span=3.3, spacing=680.0, t=12.0, mat_yield=355.0e6,
        sigma_x1=50.0, sigma_x2=50.0, sigma_y1=100.0, sigma_y2=100.0, tau_xy=5.0,
        get_puls_sp_or_up=lambda: 'UP', get_puls_boundary=lambda: 'Int', get_puls_up_boundary=lambda: up_boundary,
    )
    return SimpleNamespace(Plate=plate, Stiffener=None, E=210.0e9, v=0.3, **extra)


# ---------------------------------------------------------------------------------------------------------------------
# U3 -- AllStructure load factors are applied (pressure side recorded)
# ---------------------------------------------------------------------------------------------------------------------
def test_u3_load_factors_equal_scaling_the_inputs_by_hand():
    stf = README_BULB
    factored = flat(680, 12, 3300, stf, sx=50, sy=100, tau=5, p=0.01,
                    load_factor_stresses=1.5, load_factor_pressure=1.2).get_buckling_results()
    scaled_by_hand = flat(680, 12, 3300, stf, sx=75, sy=150, tau=7.5, p=0.012).get_buckling_results()
    unfactored = flat(680, 12, 3300, stf, sx=50, sy=100, tau=5, p=0.01).get_buckling_results()

    assert factored['valid prediction'] == scaled_by_hand['valid prediction'] == 1
    assert factored['ultimate UF raw'] == pytest.approx(scaled_by_hand['ultimate UF raw'], rel=1e-9)
    assert factored['buckling UF raw'] == pytest.approx(scaled_by_hand['buckling UF raw'], rel=1e-9)
    assert factored['ultimate UF raw'] > unfactored['ultimate UF raw']


def test_u3_load_factor_defaults_and_pressure_side_note():
    # Missing load-factor attributes mean 1.0.
    plain = S.anystructure_panel_input(_u3_structure(), lat_press=10.0)
    assert plain.axial_stress_1 == pytest.approx(50.0)
    assert plain.pressure == pytest.approx(0.010)

    structure = _u3_structure(_stress_load_factor=1.3, _lat_load_factor=2.0, _overpressure_side='stiffener side')
    factored = S.anystructure_panel_input(structure, lat_press=10.0)
    assert factored.axial_stress_1 == pytest.approx(65.0)
    assert factored.transverse_stress_2 == pytest.approx(130.0)
    assert factored.shear_stress == pytest.approx(6.5)
    assert factored.pressure == pytest.approx(0.020)

    adapter = S.solve_anystructure_panel(structure, lat_press=10.0)['diagnostics']['anystructure_adapter']
    assert adapter['stress_load_factor'] == pytest.approx(1.3)
    assert adapter['pressure_load_factor'] == pytest.approx(2.0)
    assert adapter['overpressure_side'] == 'stiffener side'
    assert any('plate side' in note for note in adapter['notes'])


# ---------------------------------------------------------------------------------------------------------------------
# U5 -- mixed-sign axial stress: both ends are solved and the envelope reported
# ---------------------------------------------------------------------------------------------------------------------
def test_u5_mixed_sign_axial_stress_reports_the_tension_end():
    mixed = flat(680, 12, 3300, README_BULB, sx=20, sx2=-300).get_buckling_results()
    tension_only = flat(680, 12, 3300, README_BULB, sx=-300).get_buckling_results()
    compression_only = flat(680, 12, 3300, README_BULB, sx=20).get_buckling_results()

    assert mixed['valid prediction'] == 1
    # A 300 MPa tension end yields at sigma / fy.
    assert mixed['ultimate UF raw'] >= 300.0 / FY - 1e-6
    assert mixed['ultimate UF raw'] == pytest.approx(
        max(tension_only['ultimate UF raw'], compression_only['ultimate UF raw']), rel=1e-12)
    assert mixed['buckling UF raw'] == pytest.approx(
        max(tension_only['buckling UF raw'], compression_only['buckling UF raw']), rel=1e-12)


def test_u5_envelope_diagnostics_name_the_governing_end():
    stiffener = SimpleNamespace(
        span=3.3, spacing=680.0, t=12.0, hw=260.0, tw=12.0, b=49.0, tf=28.0, mat_yield=355.0e6,
        sigma_x1=20.0, sigma_x2=-300.0, sigma_y1=0.0, sigma_y2=0.0, tau_xy=0.0,
        get_puls_boundary=lambda: 'Int', get_stiffener_type=lambda: 'L-bulb', get_puls_stf_end=lambda: 'C',
    )
    structure = SimpleNamespace(Plate=SimpleNamespace(get_puls_sp_or_up=lambda: 'SP'), Stiffener=stiffener,
                                E=210.0e9, v=0.3)
    labels = [label for label, _ in S.anystructure_panel_inputs(structure)]
    assert labels == ['compression-end', 'tension-end']
    result = S.solve_anystructure_panel(structure)
    envelope = result['diagnostics']['axial_stress_envelope']
    assert envelope['governing_ultimate_end'] == 'tension-end'
    assert result['ultimate_usage_factor'] == pytest.approx(300.0 / FY, rel=1e-6)
    vector = S.predict_anystructure_uf_with_acceptance(structure, selected_method='ultimate')
    assert vector[1] == pytest.approx(300.0 / FY, rel=1e-6)

    stiffener.sigma_x2 = 20.0  # same sign: one solve, no envelope
    assert [label for label, _ in S.anystructure_panel_inputs(structure)] == ['single']


@pytest.mark.parametrize('tension', [-0.5, -3.5])
def test_u5_small_tension_end_does_not_invalidate_the_panel(tension):
    # |sigma| * max_load_factor = 350 MPa < fy: the tension end cannot yield within the load range, so its
    # usage factor is below 1 / max_load_factor and the envelope is the compression end.
    mixed = flat(680, 12, 3300, README_BULB, sx=150, sx2=tension).get_buckling_results()
    compression_only = flat(680, 12, 3300, README_BULB, sx=150).get_buckling_results()
    assert abs(tension) * S.S3SolverConfig().max_load_factor < FY
    assert mixed['valid prediction'] == 1
    assert mixed['ultimate UF raw'] == pytest.approx(compression_only['ultimate UF raw'], rel=1e-12)
    assert mixed['buckling UF raw'] == pytest.approx(compression_only['buckling UF raw'], rel=1e-12)


def test_u5_envelope_buckling_diagnostics_come_from_the_buckling_end():
    # The 300 MPa tension end governs the ultimate (sigma / fy); the compression end governs the buckling UF.
    mixed = flat(900, 10, 3300, README_BULB, sx=100, sx2=-300).get_buckling_results()
    compression_only = flat(900, 10, 3300, README_BULB, sx=100).get_buckling_results()
    assert mixed['ultimate UF raw'] == pytest.approx(300.0 / FY, rel=1e-6)
    assert mixed['buckling UF raw'] == pytest.approx(compression_only['buckling UF raw'], rel=1e-12)
    assert mixed['buckling UF raw'] > mixed['ultimate UF raw']
    for key in ('elastic buckling UF raw', 'critical mode', 'critical failure family', 'controlling limit'):
        assert mixed[key] == compression_only[key]


# ---------------------------------------------------------------------------------------------------------------------
# U6 -- CSR proportions use the mm scantlings (no per-value mm/m heuristic)
# ---------------------------------------------------------------------------------------------------------------------
def test_u6_csr_check_fails_thin_8mm_plate_by_hand_calculation():
    result = flat(700, 8, 3300, (400, 8, 120, 8, 'T'), sx=50, p=0.01).get_buckling_results()
    requirement = result['CSR requirement']
    scale = math.sqrt(235.0 / FY)
    t_net = 8.0 - S.DEFAULT_ANYSTRUCTURE_CSR_CORROSION_ADDITION_MM  # 6 mm
    assert requirement['values']['plate_net_thickness'] == pytest.approx(t_net)
    # IACS CSR: s / t_net = 116.7 > 100 sqrt(235/fy) = 81.4 and h_w / t_w,net = 66.7 > 75 sqrt(235/fy) = 61.0
    assert 700.0 / t_net > 100.0 * scale
    assert requirement['checks']['plate_slenderness'] is False
    assert 400.0 / t_net > 75.0 * scale
    assert requirement['checks']['web_slenderness'] is False
    assert result['CSR'][0] == 0
    assert result['CSR color'] == 'red'


# ---------------------------------------------------------------------------------------------------------------------
# S7 -- IACS/CSR effective width coefficient 1.12 / (1 + 1.75 r^-1.6)
# ---------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('length', [1200.0, 2355.0, 3400.0, 5889.0, 9000.0])
def test_s7_effective_width_coefficient_matches_csr_formula(length):
    panel = s3_panel(length=length)
    out = S.effective_stiffener_plate_width(panel)
    ratio = (length / math.sqrt(3.0)) / 680.0
    expected = min(1.0, 1.12 / (1.0 + 1.75 * ratio ** -1.6)) if ratio >= 1.0 else 0.407 * ratio
    assert out['length_ratio'] == pytest.approx(ratio)
    assert out['coefficient'] == pytest.approx(expected, rel=1e-12)


# ---------------------------------------------------------------------------------------------------------------------
# S4 -- supports that are accepted but not modelled
# ---------------------------------------------------------------------------------------------------------------------
def test_s4_free_edge_is_rejected_and_clamped_is_flagged():
    free = S.solve_u3_panel(S.row_to_u3_input(u3_row(3000.0, 1000.0, 12.0, 100.0, rot=('SS', 'FS'))))
    assert free.valid is False
    assert free.invalid_reason == 'unsupported-rotational-support-free-edge'

    simply = S.solve_u3_panel(S.row_to_u3_input(u3_row(3000.0, 1000.0, 12.0, 100.0)))
    clamped = S.solve_u3_panel(S.row_to_u3_input(u3_row(3000.0, 1000.0, 12.0, 100.0, rot=('CL', 'CL'))))
    assert clamped.valid is True
    # CL is analysed as SS (conservative) and says so.
    assert clamped.ultimate_usage_factor == pytest.approx(simply.ultimate_usage_factor, rel=1e-12)
    assert any("'CL'" in note for note in clamped.diagnostics['support_notes'])
    assert simply.diagnostics['support_notes'] == []


def test_s4_single_free_edge_from_anystructure_is_not_turned_into_ss():
    panel = S.anystructure_panel_input(_u3_structure(up_boundary='SSSF'))
    assert 'FS' in (panel.rotational_support_1, panel.rotational_support_2)
    assert S.solve_u3_panel(panel).invalid_reason == 'unsupported-rotational-support-free-edge'


@pytest.mark.parametrize('support', ['Girder - long', 'Girder - trans'])
def test_s4_girder_in_plane_support_is_flagged_and_equals_integrated(support):
    integrated = S.solve_s3_panel(s3_panel())
    girder = S.solve_s3_panel(s3_panel(in_plane_support=support))
    assert girder.ultimate_usage_factor == pytest.approx(integrated.ultimate_usage_factor, rel=1e-12)
    assert any(support in note for note in girder.diagnostics['support_notes'])


# ---------------------------------------------------------------------------------------------------------------------
# U2 -- no lift of the ultimate capacity for pressure-dominated panels
# ---------------------------------------------------------------------------------------------------------------------
def test_u2_ultimate_uf_is_non_decreasing_in_pressure():
    stf = (200, 12, 100, 15, 'T')
    values = []
    for pressure in (0.0, 0.05, 0.10, 0.20, 0.40):
        result = flat(500, 18, 1500, stf, sx=200, p=pressure).get_buckling_results()
        assert result['valid prediction'] == 1
        values.append(result['ultimate UF raw'])
    assert all(uf >= 200.0 / FY for uf in values)
    assert all(later >= earlier - 1e-12 for earlier, later in zip(values, values[1:]))


@pytest.mark.parametrize('pressure', [0.0, 0.1, 0.4])
def test_u2_buckling_uf_is_never_below_the_ultimate_uf(pressure):
    # The buckling-strength capacity is min(ultimate, elastic buckling), so its UF bounds the ultimate UF:
    # a 'buckling' acceptance can never pass a panel whose ultimate UF fails.
    panel = flat(500, 18, 1500, (200, 12, 100, 15, 'T'), sx=200, p=pressure)
    result = panel.get_buckling_results()
    assert result['valid prediction'] == 1
    assert result['buckling UF raw'] >= result['ultimate UF raw'] - 1e-12
    panel.set_buckling_parameters(calculation_method='SemiAnalytical S3/U3', buckling_acceptance='buckling',
                                  stiffener_support='Continuous')
    assert panel.get_buckling_results()['selected UF'] >= result['ultimate UF'] - 1e-12


# ---------------------------------------------------------------------------------------------------------------------
# U1 -- ultimate capacity never above the global (column) elastic buckling load
# ---------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('stf', [(120, 10, 0, 0, 'FB'), (120, 14, 80, 12, 'T')])
def test_u1_ultimate_not_above_euler_column_with_full_plate(stf):
    result = flat(800, 16, 7000, stf, sx=100).get_buckling_results()
    hw, tw, bf, tf, _ = stf
    euler = euler_full_plate(800.0, 16.0, 7000.0, hw, tw, bf, tf)
    assert result['valid prediction'] == 1
    # The full attached plate is an upper bound of the column stiffness: UF >= sigma_x / sigma_E.
    assert result['ultimate UF raw'] >= 100.0 / euler * (1.0 - 1e-6)
    assert result['confidence'] == 'low'


def test_u1_solver_capacity_bounded_by_column_factor_and_flagged():
    panel = s3_panel(length=7000.0, stiffener_spacing=800.0, plate_thickness=16.0, stiffener_type='Flatbar',
                     stiffener_height=120.0, web_thickness=10.0, flange_width=0.0, flange_thickness=0.0,
                     axial_stress=100.0, transverse_stress_1=0.0, transverse_stress_2=0.0, shear_stress=0.0,
                     pressure=0.0)
    for config in (S.S3SolverConfig(), S.S3SolverConfig(include_solver_diagnostics=False)):
        result = S.solve_s3_panel(panel, config)
        euler = euler_full_plate(800.0, 16.0, 7000.0, 120.0, 10.0, 0.0, 0.0)
        assert result.valid is True
        assert 1.0 / result.ultimate_usage_factor <= euler / 100.0 * (1.0 + 1e-6)
        assert result.confidence == 'low'
        assert result.diagnostics['buckling_strength']['reported_ultimate_capacity_factor'] <= (
            result.diagnostics['buckling_strength']['raw_ultimate_capacity_factor'] + 1e-12)


def test_u1_ultimate_not_above_global_transverse_mode_outside_the_continuation_basis():
    # sigma_x = 0: the global orthotropic sigma_y mode (stiffeners bending over the span) is not in the
    # continuation basis, so its elastic factor bounds the ultimate (RP-C201 Sec 7 treats transverse
    # compression as stiffener bending through p0, eq. 7.9, without post-buckling reserve).
    panel = s3_panel(length=6000.0, stiffener_spacing=700.0, plate_thickness=10.0, stiffener_type='Flatbar',
                     stiffener_height=120.0, web_thickness=10.0, flange_width=0.0, flange_thickness=0.0,
                     axial_stress=0.0, transverse_stress_1=120.0, transverse_stress_2=120.0, shear_stress=0.0,
                     pressure=0.0)
    result = S.solve_s3_panel(panel)
    strength = result.diagnostics['buckling_strength']
    assert result.valid is True
    assert strength['global_column_elastic_cap']['critical_mode'] == 'stiffener-column'
    assert 1.0 / result.ultimate_usage_factor <= result.diagnostics['buckling']['stiffener_column_factor'] * (1 + 1e-12)


# ---------------------------------------------------------------------------------------------------------------------
# S1 -- the Ritz plate basis follows the aspect ratio
# ---------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('family', ['S3', 'U3'])
def test_s1_auto_basis_converges_the_ultimate(family):
    # Discretisation convergence: the auto basis (ceil(2 a/b) + 1 half-waves) is within 0.5 % of a basis with
    # ceil(3 a/b) + 1 half-waves.  ceil(a/b) + 1 half-waves was 4.5 % (S3) and 7 % (U3) low.
    if family == 'S3':
        panel, solve = s3_panel(), S.solve_s3_panel
    else:
        panel = S.U3PanelInput(length=3300.0, width=680.0, plate_thickness=12.0, yield_stress_plate=FY,
                               axial_stress_1=50.0, axial_stress_2=50.0, transverse_stress_1=100.0,
                               transverse_stress_2=100.0, shear_stress=5.0, pressure=0.01,
                               in_plane_support='Integrated', rotational_support_1='SS', rotational_support_2='SS')
        solve = S.solve_u3_panel
    fine_modes = tuple(range(1, math.ceil(3.0 * 3300.0 / 680.0) + 2))
    auto = solve(panel).ultimate_usage_factor
    fine = solve(panel, S.S3SolverConfig(longitudinal_modes=fine_modes)).ultimate_usage_factor
    assert auto == pytest.approx(fine, rel=5e-3)


@pytest.mark.parametrize('aspect', [1, 2, 3, 4, 5, 6, 8, 10, 12])
def test_s1_u3_elastic_factor_matches_classical_ss_plate(aspect):
    b, t = 1000.0, 12.0
    panel = S.row_to_u3_input(u3_row(aspect * b, b, t, 100.0))
    expected = k_uniform(aspect * b, b) * sigma_e(b, t) / 100.0
    assert u3_elastic_factor(panel) == pytest.approx(expected, rel=0.01)


@pytest.mark.parametrize('aspect', [0.2, 0.3, 0.4])
def test_s1_u3_wide_plate_under_transverse_compression(aspect):
    # sigma_y on a plate wider than long (b > a) buckles with about b/a half-waves along y.
    a, t = 1000.0, 12.0
    b = a / aspect
    panel = S.row_to_u3_input(u3_row(a, b, t, 0.0, sy1=100.0))
    expected = k_uniform(b, a) * sigma_e(a, t) / 100.0
    assert u3_elastic_factor(panel) == pytest.approx(expected, rel=0.01)


def test_s1_s3_local_plate_factor_matches_classical_ss_plate():
    # Flat bar: no local restraint factor, so the local family is the SS plate between stiffeners.
    s, t, length = 700.0, 10.0, 8.0 * 700.0 + 350.0
    panel = s3_panel(length=length, stiffener_spacing=s, plate_thickness=t, stiffener_type='Flatbar',
                     stiffener_height=200.0, web_thickness=14.0, flange_width=0.0, flange_thickness=0.0,
                     axial_stress=100.0, transverse_stress_1=0.0, transverse_stress_2=0.0, shear_stress=0.0,
                     pressure=0.0)
    config = S.S3SolverConfig()
    modes = S.build_ritz_modes(panel, S.build_section_properties(panel), config)
    local = min(mode.linear_stiffness / mode.geometric_stiffness for mode in modes
                if mode.family == 'local' and mode.geometric_stiffness > 0.0)
    assert local == pytest.approx(k_uniform(length, s) * sigma_e(s, t) / 100.0, rel=0.01)


def test_s1_s3_local_plate_wider_than_long_under_transverse_compression():
    # Stiffener spacing larger than the span: sigma_y needs about s/L half-waves across the plate.
    s, t, length = 2000.0, 14.0, 600.0
    panel = s3_panel(length=length, stiffener_spacing=s, plate_thickness=t, stiffener_type='Flatbar',
                     stiffener_height=200.0, web_thickness=14.0, flange_width=0.0, flange_thickness=0.0,
                     axial_stress=0.0, transverse_stress_1=50.0, transverse_stress_2=50.0, shear_stress=0.0,
                     pressure=0.0)
    config = S.S3SolverConfig()
    modes = S.build_ritz_modes(panel, S.build_section_properties(panel), config)
    local = min(mode.linear_stiffness / mode.geometric_stiffness for mode in modes
                if mode.family == 'local' and mode.geometric_stiffness > 0.0)
    assert local == pytest.approx(k_uniform(s, length) * sigma_e(length, t) / 50.0, rel=0.01)


def test_s1_explicit_mode_tuples_and_convergence_bases():
    panel = S.row_to_u3_input(u3_row(8000.0, 1000.0, 12.0, 100.0))
    fixed = dataclasses.replace(S.S3SolverConfig(), auto_mode_count=False)
    assert len(S.build_u3_ritz_modes(panel, fixed)) == 3 * 2
    explicit = dataclasses.replace(fixed, longitudinal_modes=tuple(range(1, 13)))
    assert u3_elastic_factor(panel, explicit) == pytest.approx(sigma_e(1000.0, 12.0) * 4.0 / 100.0, rel=0.01)

    # The mode-convergence bases must be richer than the adaptive base basis.
    checked = S.solve_u3_panel(panel, dataclasses.replace(S.S3SolverConfig(), check_mode_convergence=True))
    base = checked.diagnostics['ritz_basis']['plate_longitudinal_modes']
    medium = checked.diagnostics['mode_convergence']['medium_basis']['longitudinal_modes']
    high = checked.diagnostics['mode_convergence']['high_basis']['longitudinal_modes']
    assert max(medium) == max(base) + 1
    assert max(high) == max(base) + 2


# ---------------------------------------------------------------------------------------------------------------------
# S2 -- stress gradients drive elastic buckling (DNV-RP-C201 Sec 6.6 and 6.8)
# ---------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('aspect', [1.0, 2.0])
def test_s2_in_plane_bending_and_compression_gradient(aspect):
    b, t = 1000.0, 12.0
    reference = sigma_e(b, t)

    bending = S.row_to_u3_input(u3_row(aspect * b, b, t, 100.0, -100.0))
    k_bending = u3_elastic_factor(bending) * 100.0 / reference
    assert k_bending == pytest.approx(23.9, rel=0.05)  # classical minimum for a/b >= 2/3

    half = S.row_to_u3_input(u3_row(aspect * b, b, t, 100.0, -50.0))
    k_half = u3_elastic_factor(half) * 100.0 / reference
    assert k_half <= 13.4 * (1.0 + 1e-9)  # no unconservatism against the classical k(psi = -0.5)

    # Swapping the two edges describes the same plate.
    swapped = S.row_to_u3_input(u3_row(aspect * b, b, t, -50.0, 100.0))
    assert u3_elastic_factor(swapped) == pytest.approx(u3_elastic_factor(half), rel=1e-12)

    # Pure in-plane bending is a valid load case (was rejected as zero-variable-load);
    # the compressed edge yields no later than sigma_1 / fy.
    solved = S.solve_u3_panel(bending)
    assert solved.valid is True
    assert solved.ultimate_usage_factor >= 100.0 / FY


@pytest.mark.parametrize('aspect', [0.5, 0.75, 1.0, math.sqrt(2.0), 1.7, 2.5, 3.5])
@pytest.mark.parametrize('psi', [0.5, 0.0, -0.5, -1.0, -1.5])
def test_s2_gradient_driver_not_above_exact_ritz_solution(aspect, psi):
    b, t = 1000.0, 12.0
    panel = S.row_to_u3_input(u3_row(aspect * b, b, t, 100.0, 100.0 * psi))
    k_solver = u3_elastic_factor(panel) * 100.0 / sigma_e(b, t)
    # 1 % covers the accuracy of the RP-C201 k_sigma curve fit itself.
    assert k_solver <= k_gradient_ritz(aspect, psi) * 1.01


def test_s2_gradient_below_rp_c201_range_is_bounded_conservatively():
    # psi = -2.5 is outside RP-C201 Sec 6.6; the driver uses k_sigma(-2) = 5.98*9 (a lower bound of the
    # true k_sigma), and the tensile edge still governs yield: UF >= 250/355.
    panel = S.row_to_u3_input(u3_row(3000.0, 1000.0, 12.0, 100.0, -250.0))
    driver = S.axial_buckling_driver(panel)
    k_uniform = min((m / 3.0 + 3.0 / m) ** 2 for m in range(1, 10))
    assert driver['stress'] == pytest.approx(4 * 100.0 / (5.98 * 9) * (k_uniform / 4), rel=1e-12)
    # Yield at the tensile edge: pure tension reproduces sigma/fy exactly; with a buckling compressed part the
    # membrane redistribution lowers the edge stress by <= 0.2 % in these cases, hence the 0.5 % allowance.
    result = S.solve_u3_panel(panel)
    assert result.valid is True
    assert result.ultimate_usage_factor >= 0.995 * 250.0 / 355.0

    mixed = S.solve_u3_panel(S.row_to_u3_input(u3_row(3000.0, 1000.0, 12.0, 20.0, -300.0)))
    assert mixed.valid is True
    assert mixed.ultimate_usage_factor >= 0.995 * 300.0 / 355.0


@pytest.mark.parametrize('family', ['S3', 'U3'])
def test_s2_transverse_gradient_uses_rp_c201_6_8_design_stress(family):
    length, width, sigma_max = 3000.0, 800.0, 100.0
    # sigma_y = 100 at x = 0 and -40 at x = L.  RP-C201 Sec 6.8 by hand:
    l1 = min(0.25 * length, 0.5 * width)
    at_l1 = sigma_max + (-40.0 - sigma_max) * l1 / length
    design = max(at_l1, 0.75 * sigma_max)
    if family == 'S3':
        def build(sy1, sy2):
            return s3_panel(length=length, stiffener_spacing=width, axial_stress=0.0, shear_stress=0.0,
                            pressure=0.0, transverse_stress_1=sy1, transverse_stress_2=sy2)

        def factor(panel):
            config = S.S3SolverConfig()
            section = S.build_section_properties(panel)
            return S.elastic_buckling_factors(panel, section, S.build_ritz_modes(panel, section, config),
                                              config)['critical_factor']
    else:
        def build(sy1, sy2):
            return S.row_to_u3_input(u3_row(length, width, 12.0, 0.0, sy1=sy1, sy2=sy2))

        factor = u3_elastic_factor
    gradient = factor(build(sigma_max, -40.0))
    uniform_design = factor(build(design, design))
    assert gradient == pytest.approx(uniform_design, rel=1e-9)
    assert factor(build(-40.0, sigma_max)) == pytest.approx(gradient, rel=1e-12)


# ---------------------------------------------------------------------------------------------------------------------
# U9 -- yield checks: sigma_y varies along the length, U3 sigma_x across the width
# ---------------------------------------------------------------------------------------------------------------------
def _reference_plate_yield_ratio(panel, modes, amplitudes, field, x_values, y_values):
    total = np.asarray(amplitudes) + np.asarray([mode.imperfection for mode in modes])
    second = np.einsum('hij,i,j->h', field.coupling, total, total) - field.imperfection_offset
    gx, gy = field.grid_x_fractions, field.grid_y_fractions
    # RP-C201 Fig 6-5 / 6-6: sigma_x varies across the width (y), sigma_y along the length (x).
    sigma_x = x_values[0] + (x_values[1] - x_values[0]) * gy + field.sigma_x_grid @ second
    sigma_y = y_values[0] + (y_values[1] - y_values[0]) * gx + field.sigma_y_grid @ second
    tau = panel.shear_stress + field.tau_grid @ second
    vm = np.sqrt(np.maximum(sigma_x ** 2 - sigma_x * sigma_y + sigma_y ** 2 + 3.0 * tau ** 2, 0.0))
    return float(vm.max()) / panel.yield_stress_plate


@pytest.mark.parametrize('target', [(1, 1), (2, 2)])
def test_u9_s3_plate_yield_uses_transverse_stress_along_the_stiffener(target):
    panel = s3_panel(length=2400.0, stiffener_spacing=800.0, stiffener_type='T-bar', stiffener_height=250.0,
                     web_thickness=10.0, flange_width=100.0, flange_thickness=14.0, axial_stress=40.0,
                     transverse_stress_1=150.0, transverse_stress_2=-50.0, shear_stress=0.0, pressure=0.0)
    config = S.S3SolverConfig()
    section = S.build_section_properties(panel)
    modes = S.build_ritz_modes(panel, section, config)
    amplitudes = [6.0 if mode.family == 'local' and (mode.m, mode.n) == target else 0.0 for mode in modes]
    field = S.build_membrane_field(modes, panel.length, panel.width, panel.plate_thickness, E,
                                   [mode.imperfection for mode in modes], config.membrane_hot_spot_fractions, NU)
    expected = _reference_plate_yield_ratio(panel, modes, amplitudes, field, (40.0, 40.0), (150.0, -50.0))
    actual = S.yield_utilization(panel, section, section, modes, amplitudes, 1.0, config)['plate']
    assert actual == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize('target', [(1, 1), (2, 2)])
def test_u9_u3_plate_yield_uses_axial_stress_across_the_width(target):
    panel = S.row_to_u3_input(u3_row(2400.0, 800.0, 12.0, 150.0, -50.0, sy1=60.0, sy2=-20.0))
    config = S.S3SolverConfig()
    modes = S.build_u3_ritz_modes(panel, config)
    amplitudes = [6.0 if (mode.m, mode.n) == target else 0.0 for mode in modes]
    field = S.build_membrane_field(modes, panel.length, panel.width, panel.plate_thickness, E,
                                   [mode.imperfection for mode in modes], config.membrane_hot_spot_fractions, NU)
    expected = _reference_plate_yield_ratio(panel, modes, amplitudes, field, (150.0, -50.0), (60.0, -20.0))
    actual = S.u3_yield_utilization(panel, modes, amplitudes, 1.0, config)['plate']
    assert actual == pytest.approx(expected, rel=1e-9)
