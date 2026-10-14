'''
Regression tests for the DNV-RP-C201 (October 2010) defects listed in
docs/review/CALCULATION_REVIEW.md (findings P1, P2, P5-P16).

Expected values come from an independent implementation of the RP formulas in this
file (unstiffened plate, Sec. 6) or from section properties computed here. Where a test
needs resistance parameters of the stiffener or girder (P2, P7, P12, P16 girder), it reads
the package's intermediate values and checks the RP equation or selection built from
them; tests of input handling (P1, P8, P16) compare package runs that the RP requires to
agree or to be ordered.
'''
import copy
import math
import warnings

import pytest

import example_data as ex
from anybuckling import FlatStru
from anybuckling.prescriptive import plates as P

E, NU, FY, GAMMA_M = 210000.0, 0.3, 355.0, 1.15
GIRDER = 'Flat plate, stiffened with girder'


def _build(domain='Flat plate, stiffened', s=680, t=12, l=3300, sx=50, sx2=None, sy=100, sy2=None,
           tau=5, p=0.01, stf=(260, 12, 49, 28, 'bulb'), gir=(500, 15, 200, 25, 'T'),
           girder_span=10000, **buckling_parameters):
    panel = FlatStru(domain)
    panel.set_material(mat_yield=FY, emodule=E, material_factor=GAMMA_M, poisson=NU)
    panel.set_plate_geometry(spacing=s, thickness=t, span=l)
    panel.set_stresses(pressure=p, sigma_x1=sx, sigma_x2=sx if sx2 is None else sx2,
                       sigma_y1=sy, sigma_y2=sy if sy2 is None else sy2, tau_xy=tau)
    if domain != 'Flat plate, unstiffened':
        hw, tw, bf, tf, kind = stf
        panel.set_stiffener(hw=hw, tw=tw, bf=bf, tf=tf, stf_type=kind, spacing=s)
    if domain == GIRDER:
        hw, tw, bf, tf, kind = gir
        panel.set_girder(hw=hw, tw=tw, bf=bf, tf=tf, stf_type=kind, spacing=700, span=girder_span)
    kwargs = dict(calculation_method='DNV-RP-C201 - prescriptive', buckling_acceptance='ultimate')
    kwargs.update(buckling_parameters)
    panel.set_buckling_parameters(**kwargs)
    return panel


def _rp_unstiffened_combined(s, t, l, sx, sy, tau, psd=0.0):
    '''Independent RP-C201 Sec. 6.3-6.6 combined check for uniform stresses (compression positive).

    Longitudinal resistance by sec. 6.6, eq. (6.21)-(6.24), with psi = 1 (k_sigma = 4), the form the
    package uses for both uniform and varying stress; eq. (6.3) differs from (6.24) by 0.25 % in lambda_p.
    '''
    k_sigma = 8.2 / (1.05 + 1.0)
    lam_p = s / t / (28.4 * math.sqrt(235 / FY) * math.sqrt(k_sigma))
    cx = 1.0 if lam_p <= 0.673 else min(1.0, (lam_p - 0.055 * (3 + 1.0)) / lam_p ** 2)
    sx_rd = (cx if sx > 0 else 1.0) * FY / GAMMA_M
    lam_c = 1.1 * s / t * math.sqrt(FY / E)
    mu = 0.21 * (lam_c - 0.2)
    if lam_c <= 0.2:
        kappa = 1.0
    elif lam_c < 2.0:
        kappa = (1 + mu + lam_c ** 2 - math.sqrt((1 + mu + lam_c ** 2) ** 2 - 4 * lam_c ** 2)) / (2 * lam_c ** 2)
    else:
        kappa = 1 / (2 * lam_c ** 2) + 0.07
    h_alpha = max(0.0, 0.05 * s / t - 0.75)
    kp = 1.0 if psd <= 2 * (t / s) ** 2 * FY else max(0.0, 1.0 - h_alpha * (psd / FY - 2 * (t / s) ** 2))
    a = 1.3 * t / l * math.sqrt(E / FY)
    sy_rd = ((a + kappa * (1 - a)) * FY * kp if sy > 0 else FY) / GAMMA_M
    kl = 5.34 + 4 * (s / l) ** 2 if l >= s else 5.34 * (s / l) ** 2 + 4
    lam_w = 0.795 * s / t * math.sqrt(FY / (E * kl))
    if sy > 0:      # (6.19)/(6.20): elastic shear resistance when sigma_y is compressive
        c_tau = 1.0 if lam_w <= 0.8 else (1 - 0.8 * (lam_w - 0.8) if lam_w <= 1.25 else 1 / lam_w ** 2)
    else:           # (6.14)/(6.15)
        c_tau = 1.0 if lam_w <= 0.8 else (1 - 0.625 * (lam_w - 0.8) if lam_w <= 1.2 else 0.9 / lam_w)
    tau_rd = c_tau * FY / (GAMMA_M * math.sqrt(3))
    if sx > 0 and sy > 0:
        ci = 1 - s / (120 * t) if s / t <= 120 else 0.0
    else:
        ci = 1.0
    rx, ry = sx / sx_rd, sy / sy_rd
    return {'combined': rx ** 2 + ry ** 2 - ci * rx * ry + (tau / tau_rd) ** 2, 'lambda_w': lam_w,
            'sy_rd': sy_rd, 'tau_rd': tau_rd}


def _composite(parts):
    '''Area, centroid and second moment of a list of rectangles (width, height, z_centre).'''
    area = sum(b * h for b, h, _ in parts)
    zc = sum(b * h * z for b, h, z in parts) / area
    inertia = sum(b * h ** 3 / 12 + b * h * (z - zc) ** 2 for b, h, z in parts)
    return area, zc, inertia


# --------------------------------------------------------------------------------------------
# P1 - girder span must be settable through the public API
# --------------------------------------------------------------------------------------------
def test_girder_span_is_settable_and_used():
    short = _build(GIRDER, sy=40, girder_span=4000).get_buckling_results()['Girder']
    long_ = _build(GIRDER, sy=40, girder_span=16000).get_buckling_results()['Girder']

    reference = _build(GIRDER, sy=40, girder_span=10000)
    reference._FlatStructure.Plate.girder_lg = 4.0
    reference._FlatStructure.Stiffener.girder_lg = 4.0
    reference._FlatStructure.Girder.girder_lg = 4.0

    assert short == pytest.approx(reference.get_buckling_results()['Girder'])
    assert long_['Overpressure plate side'] > 2 * short['Overpressure plate side']


def test_girder_without_span_warns_about_placeholder():
    panel = FlatStru(GIRDER)
    panel.set_material(mat_yield=FY, emodule=E, material_factor=GAMMA_M, poisson=NU)
    panel.set_plate_geometry(spacing=680, thickness=12, span=3300)
    panel.set_stresses(pressure=0.01, sigma_x1=50, sigma_x2=50, sigma_y1=40, sigma_y2=40, tau_xy=5)
    with pytest.warns(UserWarning, match='girder span'):
        panel.set_girder(hw=500, tw=15, bf=200, tf=25, stf_type='T', spacing=700)


# --------------------------------------------------------------------------------------------
# P2 - sniped girders: (7.63)/(7.64) apply when q*LG^2/8 < N*z*
# --------------------------------------------------------------------------------------------
def test_sniped_girder_uses_axial_dominated_equations():
    structure = _build(GIRDER, sx=0, sy=40, tau=0, p=0.0, girder_support='Sniped')._FlatStructure
    unstiffened = structure.unstiffened_plate_buckling()
    stiffened = structure.stiffened_panel(unstf_pl_data=unstiffened)
    data = structure.girder(unstf_pl_data=unstiffened, stf_pl_data=stiffened)

    n, z, q, lg = data['NySd'], data['zstar'], data['qSd girder side'], data['LG']
    amplification = 1 - n / data['NE']
    assert q * lg ** 2 / 8 < n * z
    eq_7_63 = n / data['NksRd'] - 2 * n / data['NRd'] + (n * z - q * lg ** 2 / 8) / (data['MstRd'] * amplification)
    eq_7_64 = n / data['NkpRd'] + (n * z - q * lg ** 2 / 8) / (data['MpRd'] * amplification)
    assert data['UF Simplified girder side'] == pytest.approx(max(0.0, eq_7_63, eq_7_64))
    assert data['UF Simplified girder side'] > 1.0


# --------------------------------------------------------------------------------------------
# P5, P9, P14, P15 - unstiffened plate checks (Sec. 6)
# --------------------------------------------------------------------------------------------
@pytest.mark.parametrize('s, t, l, sx, sy, tau', [
    (1000, 10, 3000, 30, 20, 40),     # P5: sigma_y compressive, lambda_w > 1.25 -> C_tau_e
    (800, 10, 2400, 60, 40, 60),      # P5: 0.8 < lambda_w <= 1.25
    (800, 10, 2400, 60, -40, 60),     # P9 + P14: sigma_y tensile -> C_tau (0.625 form), ci = 1
    (800, 9, 2400, -60, 40, 20),      # P14: sigma_x tensile -> ci = 1
])
def test_unstiffened_combined_check_matches_rp(s, t, l, sx, sy, tau):
    panel = _build('Flat plate, unstiffened', s=s, t=t, l=l, sx=sx, sy=sy, tau=tau, p=0.0)
    code = panel._FlatStructure.unstiffened_plate_buckling()
    reference = _rp_unstiffened_combined(s, t, l, sx, sy, tau)

    assert code['UF Combined stresses'] == pytest.approx(reference['combined'], rel=1e-9)


def test_shear_reduction_factor_is_continuous_and_matches_rp():
    # lambda_w = 1.2 boundary: (6.15) gives 0.75 from both sides
    for s, t, l in ((800, 9, 2400), (700, 9, 2100), (900, 9, 2700)):
        panel = _build('Flat plate, unstiffened', s=s, t=t, l=l, sx=0, sy=0, tau=40, p=0.0)
        code = panel._FlatStructure.unstiffened_plate_buckling()
        reference = _rp_unstiffened_combined(s, t, l, 0.0, 0.0, 40.0)
        assert code['UF Shear stresses'] == pytest.approx(40.0 / reference['tau_rd'], rel=1e-9)


def test_lateral_pressure_reduction_factor_is_not_negative():
    # (6.10): kp >= 0; a negative kp would give a negative resistance and a negative (passing) UF
    # h_alpha = 0.05*125 - 0.75 = 5.5 -> unclamped kp = 1 - 5.5*(100/355 - 2*(8/1000)**2) = -0.55
    panel = _build('Flat plate, unstiffened', s=1000, t=8, l=3000, sx=0, sy=50, tau=0, p=100.0)
    code = panel._FlatStructure.unstiffened_plate_buckling()

    assert code['syR'] == 0
    # zero resistance against a compressive transverse stress is a failure, not a pass
    assert code['UF transverse stresses'] == math.inf
    assert code['UF Combined stresses'] == math.inf


# --------------------------------------------------------------------------------------------
# P6 - girder section moduli use le <= 0.3 LG (continuous)
# --------------------------------------------------------------------------------------------
def test_girder_section_moduli_respect_shear_lag_limit():
    structure = _build(GIRDER, sx=10, sy=10, tau=5, p=0.002, girder_span=6000)._FlatStructure
    unstiffened = structure.unstiffened_plate_buckling()
    stiffened = structure.stiffened_panel(unstf_pl_data=unstiffened)
    data = structure.girder(unstf_pl_data=unstiffened, stf_pl_data=stiffened)

    assert data['le'] > 0.3 * 6000
    le = 0.3 * 6000
    t, hw, tw, bf, tf = 12.0, 500.0, 15.0, 200.0, 25.0
    _, zc, inertia = _composite([(le, t, 0.0), (tw, hw, t / 2 + hw / 2), (bf, tf, t / 2 + hw + tf / 2)])
    assert data['Wep'] == pytest.approx(inertia / zc, rel=1e-9)
    assert data['WeG'] == pytest.approx(inertia / (t / 2 + hw + tf - zc), rel=1e-9)


# --------------------------------------------------------------------------------------------
# P7 - p0 for tensile sigma_x, eq. (8.7)
# --------------------------------------------------------------------------------------------
def test_girder_p0_for_tensile_longitudinal_stress():
    structure = _build(GIRDER, sx=-60, sy=40, tau=80, panel_length_Lp=6 * 3300)._FlatStructure
    unstiffened = structure.unstiffened_plate_buckling()
    stiffened = structure.stiffened_panel(unstf_pl_data=unstiffened)
    data = structure.girder(unstf_pl_data=unstiffened, stf_pl_data=stiffened)

    t, s, l, lg, hw = 12.0, 680.0, 3300.0, 10000.0, 500.0
    area_stiffener = 260 * 12 + 49 * 28
    assert data['C'] > 0
    expected = 0.4 * (t + area_stiffener / s) / (hw * (1 - s / lg)) * FY / E * (lg / l) ** 2 * data['C'] * 80
    assert data['p0'] == pytest.approx(expected, rel=1e-9)


# --------------------------------------------------------------------------------------------
# P8 - unknown panel length must not give the most favourable tau_crg
# --------------------------------------------------------------------------------------------
def test_unknown_panel_length_is_not_less_conservative_than_any_finite_length():
    unknown = _build(GIRDER, sx=60, sy=40, tau=80).get_buckling_results()['Girder']
    for spans in (3, 6, 12, 50):
        finite = _build(GIRDER, sx=60, sy=40, tau=80, panel_length_Lp=spans * 3300).get_buckling_results()['Girder']
        assert unknown['Overpressure plate side'] >= finite['Overpressure plate side'] - 1e-12


# --------------------------------------------------------------------------------------------
# P10 - V_Sd >= V_Rd must report failure, not raise
# --------------------------------------------------------------------------------------------
def test_stiffener_shear_overload_reports_failure_instead_of_raising():
    results = _build(s=680, t=12, l=3300, sx=50, sy=20, tau=5, p=0.20,
                     stf=(200, 6, 0, 0, 'FB')).get_buckling_results()['Stiffener']

    assert results['Shear capacity'] > 1
    assert results['Overpressure plate side'] > 1
    assert results['Overpressure stiffener side'] > 1


# --------------------------------------------------------------------------------------------
# P11 - special provisions: transverse stress sampled at l1 from the most stressed end
# --------------------------------------------------------------------------------------------
def test_special_provisions_do_not_depend_on_end_labelling():
    data = copy.deepcopy(ex.obj_dict)
    data['panel or shell'] = ['panel', '']
    scantling = P.CalcScantlings(data)
    scantling.set_stresses(sigy1=200, sigy2=0, sigx1=50, sigx2=50, tauxy=0)
    first = (scantling.get_dnv_min_section_modulus(200), scantling.get_dnv_min_thickness(200))
    scantling.set_stresses(sigy1=0, sigy2=200, sigx1=50, sigx2=50, tauxy=0)
    second = (scantling.get_dnv_min_section_modulus(200), scantling.get_dnv_min_thickness(200))

    assert first == pytest.approx(second)
    # stress at l1 = min(0.25*3.3, 0.5*0.68) = 0.34 m from the 200 MPa end: 200*(1-0.34/3.3)
    sigma_y = 200 * (1 - 0.34 / 3.3)
    sigma_j = math.sqrt(50 ** 2 + sigma_y ** 2 - 50 * sigma_y)
    fyd = 355 / 1.10
    z_expected = 3.3 ** 2 * 0.68 * 200 / (12 * (fyd - sigma_j) * 1) * 1e6 / 1e9
    assert first[0] == pytest.approx(z_expected, rel=1e-9)


# --------------------------------------------------------------------------------------------
# P12 - p0 uses psi = sigma_y2 / sigma_y1, eq. (7.9)
# --------------------------------------------------------------------------------------------
def test_equivalent_lateral_pressure_uses_transverse_stress_ratio():
    structure = _build(sx=50, sy=80, sy2=40, tau=5, p=0.01)._FlatStructure
    unstiffened = structure.unstiffened_plate_buckling()
    data = structure.stiffened_panel(unstf_pl_data=unstiffened)

    assert data['p0'] == pytest.approx((0.6 + 0.4 * 0.5) * data['C0'] * 80, rel=1e-12)


# --------------------------------------------------------------------------------------------
# P13 - plastic section modulus with the plastic neutral axis in the plate
# --------------------------------------------------------------------------------------------
@pytest.mark.parametrize('spacing, plate', [(0.68, 0.025), (0.15, 0.008)])
def test_plastic_section_modulus_matches_independent_calculation(spacing, plate):
    data = copy.deepcopy(ex.obj_dict)
    data['panel or shell'] = ['panel', '']
    data['spacing'] = [spacing, 'm']
    data['plate_thk'] = [plate, 'm']
    structure = P.Structure(data)
    hw, tw, bf, tf = 0.250297358, 0.012, 0.052, 0.029702642
    parts = [(spacing, plate, plate / 2), (tw, hw, plate + hw / 2), (bf, tf, plate + hw + tf / 2)]

    area = sum(b * h for b, h, _ in parts)
    lo, hi = 0.0, plate + hw + tf
    for _ in range(200):     # bisection for the equal-area axis
        z = 0.5 * (lo + hi)
        below = sum(b * min(max(z - (zc - h / 2), 0.0), h) for b, h, zc in parts)
        lo, hi = (z, hi) if below < area / 2 else (lo, z)
    z = 0.5 * (lo + hi)
    expected = 0.0
    for b, h, zc in parts:
        bottom, top = zc - h / 2, zc + h / 2
        if top <= z or bottom >= z:
            expected += b * h * abs(zc - z)
        else:
            expected += b * (z - bottom) ** 2 / 2 + b * (top - z) ** 2 / 2

    assert structure.get_plasic_section_modulus() == pytest.approx(expected, rel=1e-9)


# --------------------------------------------------------------------------------------------
# P16 - l_T for M_s1,Rd / M_s2,Rd is the smaller of 0.4 l / 0.8 l and the lateral support distance
# --------------------------------------------------------------------------------------------
def test_lateral_support_distance_equal_to_span_matches_no_lateral_support():
    without = _build(sx=150, sy=30, p=0.05, stf=(400, 12, 150, 20, 'T')).get_buckling_results()['Stiffener']
    with_span = _build(sx=150, sy=30, p=0.05, stf=(400, 12, 150, 20, 'T'),
                       stf_dist_between_lateral_supp=3300).get_buckling_results()['Stiffener']

    assert with_span['Overpressure plate side'] == pytest.approx(without['Overpressure plate side'], rel=1e-9)
    assert with_span['Overpressure stiffener side'] == pytest.approx(without['Overpressure stiffener side'],
                                                                     rel=1e-9)


def test_zero_transverse_resistance_does_not_crash_stiffened_panel():
    results = _build(s=1000, t=8, l=3000, sx=50, sy=50, tau=0, p=70.0,
                     stf=(300, 12, 100, 15, 'T')).get_buckling_results()

    assert results['Stiffener']['Resistance between stiffeners'] == math.inf


def test_girder_span_is_kept_on_repeated_set_girder():
    panel = _build(GIRDER, sy=40, girder_span=4000)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        panel.set_girder(hw=600, tw=15, bf=200, tf=25, stf_type='T', spacing=700)

    assert panel._FlatStructure.Plate.girder_lg == pytest.approx(4.0)


def test_girder_lateral_support_closer_than_0_4_lg_governs_both_moment_resistances():
    # eq. (8.12)/(8.13): l_T = min(0.4 L_G, L_GT) and min(0.8 L_G, L_GT); with L_GT = 2 m < 0.4 L_G
    # both use 2 m, so M_s1,Rd = M_s2,Rd
    structure = _build(GIRDER, sx=30, sy=60, tau=5, p=0.02, gir=(500, 15, 100, 15, 'L'),
                       girder_dist_between_lateral_supp=2000)._FlatStructure
    unstiffened = structure.unstiffened_plate_buckling()
    stiffened = structure.stiffened_panel(unstf_pl_data=unstiffened)
    data = structure.girder(unstf_pl_data=unstiffened, stf_pl_data=stiffened)

    assert data['Ms1Rd'] == pytest.approx(data['Ms2Rd'], rel=1e-12)


def test_girder_shear_overload_reports_failure():
    results = _build(GIRDER, sx=30, sy=40, tau=5, p=0.03,
                     gir=(500, 4, 200, 25, 'T')).get_buckling_results()['Girder']

    assert results['Shear capacity'] > 1
    assert results['Overpressure plate side'] == math.inf
    assert results['Overpressure girder side'] == math.inf


def test_no_equivalent_pressure_when_transverse_stress_is_nowhere_compressive():
    # sigma_y = 0 at one end and tensile at the other: sigma_y1,Sd = 0 in eq. (7.9) -> p0 = 0
    structure = _build(sx=80, sy=0, sy2=-100, tau=0, p=0.01, stf=(400, 12, 150, 20, 'T'))._FlatStructure
    unstiffened = structure.unstiffened_plate_buckling()
    data = structure.stiffened_panel(unstf_pl_data=unstiffened)

    assert data['p0'] == 0
