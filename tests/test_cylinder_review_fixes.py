'''
Regression tests for the DNV-RP-C202 (January 2013) review fixes C1-C10, C12-C14 and C18
(docs/review/CALCULATION_REVIEW.md, section 4.2).

Every expected value comes from an independent hand calculation written here from the equations of
the RP, or from an invariance property the standard requires (sign of the bending stress for a full
cylinder, continuity, a degenerate cone equal to a cylinder). None of them is taken from the code.
'''
import math

import pytest

from anybuckling import CylStru
from anybuckling.conversions import helper_cylinder_stress_to_force_to_stress
from anybuckling.prescriptive.cylinders import CylinderAndCurvedPlate

E, NU, FY = 210000.0, 0.3, 355.0
D = math.pi ** 2 * E / (12 * (1 - NU ** 2))  # pi^2 E / (12 (1 - nu^2)) of (3.3.1), (3.4.1), (3.6.3)

LONG = (232, 12, 49, 28, 'bulb', 680.367)
RING = (260, 12, 49, 28, 'bulb', 680)
GIRDER = (500, 15, 200, 25, 'T', 3300)


def build(domain='Unstiffened shell', r=6500, t=19, l=3300, L=20000, long=None, ring=None, girder=None,
          panel=680, mf=1.15, fab=('Fabricated', 'Fabricated'), endcap=True, k=1.0):
    c = CylStru(calculation_domain=domain)
    c.set_material(mat_yield=FY, emodule=E, material_factor=mf, poisson=NU)
    c.set_imperfection()
    c.set_fabrication_method(*fab)
    c.set_end_cap_pressure_included_in_stress(endcap)
    c.set_uls_or_als()
    c.set_length_between_girder(val=3300)
    c.set_panel_spacing(val=panel)
    c.set_shell_geometry(radius=r, thickness=t, tot_length_of_shell=L, distance_between_rings=l)
    if long:
        c.set_longitudinal_stiffener(*long)
    if ring:
        c.set_ring_stiffener(*ring)
    if girder:
        c.set_exclude_ring_stiffener()
        c.set_ring_girder(*girder)
    c.set_shell_buckling_parmeters(k)
    return c


def build_cone(r1=4000, r2=6500, length=5000, t=20, endcap=True):
    c = CylStru('Unstiffened conical shell')
    c.set_material(mat_yield=FY, emodule=E, material_factor=1.15, poisson=NU)
    c.set_imperfection()
    c.set_fabrication_method()
    c.set_end_cap_pressure_included_in_stress(endcap)
    c.set_uls_or_als()
    c.set_conical_shell_geometry(r1=r1, r2=r2, length=length, thickness=t)
    return c


# ---------------------------------------------------------------------------------------------------------
# Hand calculations straight from DNV-RP-C202 (2013)
# ---------------------------------------------------------------------------------------------------------
def reduced_coefficient(psi, xi, rho):
    return psi * math.sqrt(1 + (rho * xi / psi) ** 2)                          # (3.3.2), (3.4.2)


def gamma_m(lam, mf=1.15):
    '''(3.1.3); for other material factors the whole 1.15-based curve is scaled by mf/1.15.'''
    if lam < 0.5:
        return mf
    return (0.85 + 0.6 * lam if lam <= 1.0 else 1.45) * mf / 1.15


def fks(lam):
    return FY / math.sqrt(1 + lam ** 4)                                         # (3.2.1)


def cylinder_elastic_strengths(r, t, l, end_cap_in_axial_stress=True):
    '''Unstiffened circular cylinder, 3.4.2, Table 3-2 and (3.4.4)/(3.4.5).'''
    Zl = l ** 2 / (r * t) * math.sqrt(1 - NU ** 2)                              # (3.4.3)
    k = D * (t / l) ** 2                                                        # (3.4.1)
    fEa = reduced_coefficient(1, 0.702 * Zl, 0.5 * (1 + r / (150 * t)) ** -0.5) * k
    fEm = reduced_coefficient(1, 0.702 * Zl, 0.5 * (1 + r / (300 * t)) ** -0.5) * k
    if l / r > 3.85 * math.sqrt(r / t):
        fEt = 0.25 * E * (t / r) ** 1.5                                         # (3.4.4)
    else:
        fEt = reduced_coefficient(5.34, 0.856 * Zl ** 0.75, 0.6) * k
    if l / r > 2.25 * math.sqrt(r / t):
        fEh = 0.25 * E * (t / r) ** 2                                           # (3.4.5)
    else:
        # Table 3-2 note 1 (lateral, psi = 4) / note 2 (hydrostatic, psi = 2: axial stress excludes end cap)
        fEh = reduced_coefficient(4 if end_cap_in_axial_stress else 2, 1.04 * math.sqrt(Zl), 0.6) * k
    return fEa, fEm, fEt, fEh


def hand_cylinder_shell(r, t, l, sa=0.0, sm=0.0, tau=0.0, sh=0.0, end_cap_in_axial_stress=True, mf=1.15):
    '''Shell buckling of an unstiffened circular cylinder (3.2, 3.4); returns UF, lambda_s, gammaM.'''
    fEa, fEm, fEt, fEh = cylinder_elastic_strengths(r, t, l, end_cap_in_axial_stress)
    sx = min(sa, sa + sm, sa - sm)                                              # most compressed generator
    sj = math.sqrt(sx ** 2 - sx * sh + sh ** 2 + 3 * tau ** 2)                 # (3.2.3)
    lam2 = FY / sj * (max(0.0, -sa) / fEa + abs(sm) / fEm + max(0.0, -sh) / fEh + tau / fEt)  # (3.2.2)
    lam = math.sqrt(lam2)
    gm = gamma_m(lam, mf)
    return sj / (fks(lam) / gm), lam, gm


def hand_curved_panel_uf(r, t, l, s, sa=0.0, tau=0.0, sh=0.0):
    '''Unstiffened curved panel (3.3.2, Table 3-1), panel width s < l, no bending.'''
    Zs = s ** 2 / (r * t) * math.sqrt(1 - NU ** 2)                              # (3.3.3)
    k = D * (t / s) ** 2                                                        # (3.3.1)
    fEa = reduced_coefficient(4, 0.702 * Zs, 0.5 * (1 + r / (150 * t)) ** -0.5) * k
    fEt = reduced_coefficient(5.34 + 4 * (s / l) ** 2, 0.856 * math.sqrt(s / l) * Zs ** 0.75, 0.6) * k
    fEh = reduced_coefficient((1 + (s / l) ** 2) ** 2, 1.04 * (s / l) * math.sqrt(Zs), 0.6) * k
    sj = math.sqrt(sa ** 2 - sa * sh + sh ** 2 + 3 * tau ** 2)
    lam = math.sqrt(FY / sj * (max(0.0, -sa) / fEa + max(0.0, -sh) / fEh + tau / fEt))
    return sj / (fks(lam) / gamma_m(lam))


def ring_parameters(r, t, l, area_ring):
    '''(2.2.10)-(2.2.13): alpha, zeta and l_eo of a ring-stiffened cylinder.'''
    beta = l / (1.56 * math.sqrt(r * t))
    leo = l / beta * (math.cosh(2 * beta) - math.cos(2 * beta)) / (math.sinh(2 * beta) + math.sin(2 * beta))
    zeta = max(0.0, 2 * (math.sinh(beta) * math.cos(beta) + math.cosh(beta) * math.sin(beta)) /
               (math.sinh(2 * beta) + math.sin(2 * beta)))
    return area_ring / (leo * t), zeta, leo


def column_curve(fak, lam):
    return (1 - 0.28 * lam ** 2) * fak if lam <= 1.34 else 0.9 * fak / lam ** 2    # (3.8.5), (3.8.6)


def hand_unstiffened_column_uf(r, t, l, Lc, sa=0.0, sm=0.0, sh=0.0, k=1.0):
    '''Column buckling of a bare tube (3.8) with fak from (3.8.8)-(3.8.11); gammaM of the shell check.'''
    fEa, _, _, fEh = cylinder_elastic_strengths(r, t, l)
    a = 1 + FY ** 2 / fEa ** 2                                                  # (3.8.9)
    b = (2 * FY ** 2 / (fEa * fEh) - 1) * sh                                    # (3.8.10)
    c = sh ** 2 + FY ** 2 * sh ** 2 / fEh ** 2 - FY ** 2                        # (3.8.11)
    fak = (b + math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)                         # (3.8.8)
    Ac, Ic = 2 * math.pi * r * t, math.pi * r ** 3 * t
    ic = math.sqrt(Ic / Ac)
    fE = math.pi ** 2 * E * Ic / ((k * Lc) ** 2 * Ac)                           # (3.8.3)
    lam = math.sqrt(fak / fE)                                                   # (3.8.7)
    _, _, gm = hand_cylinder_shell(r, t, l, sa=sa, sm=sm, sh=sh)
    sa0 = max(0.0, -sa)
    uf = sa0 / (column_curve(fak, lam) / gm) + (abs(sm) / (1 - sa0 / fE)) / (fak / gm)  # (3.8.2)
    need = (k * Lc / ic) ** 2 >= 2.5 * E / FY                                   # (3.8.1)
    return uf, lam, need


# ---------------------------------------------------------------------------------------------------------
# C1 - additional hoop stress for ring-stiffened shells
# ---------------------------------------------------------------------------------------------------------
def test_c1_ring_stiffened_shell_uses_additional_hoop_stress():
    # With p = 0 and sigma_x = 0, (2.2.9) leaves only the additional hoop stress, so the shell between the
    # rings must give exactly the unstiffened-cylinder result for the same stresses.
    # (endcap=False keeps this test independent of fix C13: psi = 2 before and after.)
    ring = build('Ring Stiffened shell', ring=RING, endcap=False)
    ring.set_stresses(sasd=0, tQsd=10, shsd=-50)
    plain = build('Unstiffened shell', endcap=False)
    plain.set_stresses(sasd=0, tQsd=10, shsd=-50)

    hand_uf, _, _ = hand_cylinder_shell(6500, 19, 3300, tau=10, sh=-50, end_cap_in_axial_stress=False)
    assert plain.get_buckling_results()['Unstiffened shell'] == pytest.approx(hand_uf, rel=1e-9)
    assert ring.get_buckling_results()['Unstiffened shell'] == pytest.approx(hand_uf, rel=1e-9)


def test_c1_ring_stiffened_hoop_stress_matches_eq_2_2_9_plus_additional():
    c = build('Ring Stiffened shell', ring=RING)
    c.set_stresses(sasd=-40, tQsd=10, psd=-0.1, shsd=-50)
    r, t, l = 6500.0, 19.0, 3300.0
    alpha, zeta, _ = ring_parameters(r, t, l, area_ring=260 * 12 + 49 * 28)
    p, sx = -0.1, -40.0
    expected = p * r / t - alpha * zeta / (alpha + 1) * (p * r / t - NU * sx) - 50  # (2.2.9) + additional
    assert c._CylinderMain.shell_buckling()['shsd'][1] == pytest.approx(expected, rel=1e-9)


def test_c1_longitudinally_stiffened_hoop_stress_not_double_counted():
    c = build('Longitudinal Stiffened shell', long=LONG)
    c.set_stresses(sasd=-40, psd=-0.1, shsd=-50)
    expected = -0.1 * 6500 / 19 - 50                                           # (2.2.8) + additional
    shsd = c._CylinderMain.shell_buckling()['shsd']
    assert shsd[0] == pytest.approx(expected, rel=1e-12)
    assert shsd[1] == pytest.approx(expected, rel=1e-12)


# ---------------------------------------------------------------------------------------------------------
# C2 - conical shell: stress-input hoop stress
# ---------------------------------------------------------------------------------------------------------
def test_c2_cone_stress_input_keeps_additional_hoop_stress():
    # A cone with r1 = r2 (alpha = 0) is the cylinder r = r1, l = length (4.3.1)/(4.3.2), te = t.
    # (endcap=False keeps the C2 tests independent of fix C13: psi = 2 before and after.)
    cone = build_cone(r1=4000, r2=4000, length=5000, t=20, endcap=False)
    cone.set_conical_stresses(sasd=-60, shsd=-50)
    uf_cone = cone.get_buckling_results()['Unstiffened conical shell']

    hand_uf, _, _ = hand_cylinder_shell(4000, 20, 5000, sa=-60, sh=-50, end_cap_in_axial_stress=False)
    assert uf_cone == pytest.approx(hand_uf, rel=1e-9)

    no_hoop = build_cone(r1=4000, r2=4000, length=5000, t=20, endcap=False)
    no_hoop.set_conical_stresses(sasd=-60, shsd=0)
    assert no_hoop.get_buckling_results()['Unstiffened conical shell'] == pytest.approx(
        hand_cylinder_shell(4000, 20, 5000, sa=-60, end_cap_in_axial_stress=False)[0], rel=1e-9)


def test_c2_cone_force_input_hoop_stress_not_double_counted():
    # Force input: sigma_a = p r/(2 te) + N/(2 pi r te) (4.2.2), sigma_h = p r/te (4.2.4), nothing else.
    r, t, p, N = 4000.0, 20.0, -0.2, -2000.0
    cone = build_cone(r1=r, r2=r, length=5000, t=t, endcap=False)
    cone.set_conical_forces(Nsd=N, psd=p)
    sa = p * r / (2 * t) + N * 1000 / (2 * math.pi * r * t)
    hand_uf, _, _ = hand_cylinder_shell(r, t, 5000, sa=sa, sh=p * r / t, end_cap_in_axial_stress=False)
    assert cone.get_buckling_results()['Unstiffened conical shell'] == pytest.approx(hand_uf, rel=1e-9)

    # Stress input with the same pressure and no additional hoop stress gives the same result.
    stress = build_cone(r1=r, r2=r, length=5000, t=t, endcap=False)
    stress.set_conical_stresses(sasd=sa, psd=p, shsd=0)
    assert stress.get_buckling_results()['Unstiffened conical shell'] == pytest.approx(hand_uf, rel=1e-9)


def test_c2_cone_force_input_survives_later_geometry_change():
    changed = build_cone(t=20)
    changed.set_conical_forces(Nsd=-2000, M1sd=1000, psd=-0.2)
    changed.set_conical_shell_geometry(r1=4000, r2=6500, length=5000, thickness=25)
    fresh = build_cone(t=25)
    fresh.set_conical_forces(Nsd=-2000, M1sd=1000, psd=-0.2)
    assert changed.get_buckling_results()['Unstiffened conical shell'] == pytest.approx(
        fresh.get_buckling_results()['Unstiffened conical shell'], rel=1e-12)


# ---------------------------------------------------------------------------------------------------------
# C3 - longitudinally stiffened shell at sigma_x = 0
# ---------------------------------------------------------------------------------------------------------
def test_c3_longitudinally_stiffened_shell_continuous_at_zero_axial_stress():
    def uf(sa):
        c = build('Longitudinal Stiffened shell', long=LONG)
        c.set_stresses(sasd=sa, tQsd=50, shsd=-100)
        return c.get_buckling_results()['Longitudinal stiffened shell']

    # se -> 0 when sigma_x -> 0 (3.6.7), C -> rho*xi (3.6.4): the hoop and shear terms of (3.2.2) remain.
    assert uf(0) == pytest.approx(uf(-1e-4), rel=1e-4)
    assert uf(0) == pytest.approx(uf(-0.5), rel=5e-3)


# ---------------------------------------------------------------------------------------------------------
# C4 - full longitudinally stiffened shells see +/- sigma_m
# ---------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('domain, girder', [('Longitudinal Stiffened shell', None),
                                            ('Orthogonally Stiffened shell', GIRDER)])
def test_c4_full_longitudinally_stiffened_shell_independent_of_bending_sign(domain, girder):
    def uf(sm):
        c = build(domain, long=LONG, girder=girder)
        c.set_stresses(sasd=-20, smsd=sm, shsd=-100)
        return c.get_buckling_results()['Longitudinal stiffened shell']

    assert uf(50) == pytest.approx(uf(-50), rel=1e-12)


# ---------------------------------------------------------------------------------------------------------
# C5 / C6 - column buckling
# ---------------------------------------------------------------------------------------------------------
def test_c5_unstiffened_tube_column_buckling_axial():
    c = build('Unstiffened shell', r=500, t=12, l=3000, L=40000)
    c.set_stresses(sasd=-150)
    res = c.get_buckling_results()

    hand_uf, lam, need = hand_unstiffened_column_uf(500, 12, 3000, 40000, sa=-150)
    assert lam > 1.34 and need
    assert res['Need to check column buckling'] is True
    assert res['Column stability UF'] == pytest.approx(hand_uf, rel=1e-9)
    assert hand_uf == pytest.approx(1.1917, abs=1e-4)
    assert res['Column stability check'] is False


def test_c5_unstiffened_tube_column_buckling_with_bending_and_pressure():
    p = -0.3
    c = build('Unstiffened shell', r=500, t=12, l=3000, L=40000)
    c.set_stresses(sasd=-60, smsd=40, psd=p)
    hand_uf, lam, _ = hand_unstiffened_column_uf(500, 12, 3000, 40000, sa=-60, sm=40, sh=p * 500 / 12)
    assert lam > 1.34
    assert c.get_buckling_results()['Column stability UF'] == pytest.approx(hand_uf, rel=1e-9)


def test_c5_unstiffened_panel_reports_column_buckling():
    # Panel wider than l: reported as a circular cylinder (3.3.2), so the column check uses its gammaM.
    c = build('Unstiffened panel', r=500, t=12, l=3000, L=40000, panel=4000)
    c.set_stresses(sasd=-150)
    hand_uf, _, _ = hand_unstiffened_column_uf(500, 12, 3000, 40000, sa=-150)
    assert c.get_buckling_results()['Column stability UF'] == pytest.approx(hand_uf, rel=1e-9)


def test_c6_column_curve_beyond_lambda_1_34():
    f = CylinderAndCurvedPlate._column_reduced_buckling_strength
    for lam in (1.3401, 1.5, 2.0, 3.0):
        assert f(1.0, lam) == pytest.approx(0.9 / lam ** 2, rel=1e-12)
    # (3.8.5) and (3.8.6) meet within 1 % of fak at lambda = 1.34
    assert abs(f(1.0, 1.34 + 1e-9) - f(1.0, 1.34)) < 0.01


# ---------------------------------------------------------------------------------------------------------
# C7 - force/stress conversion units
# ---------------------------------------------------------------------------------------------------------
def test_c7_set_forces_gives_membrane_stresses_in_pa():
    N, M, T, Q = -1000.0, 10000.0, 10000.0, 500.0                               # kN, kNm, kNm, kN
    c = build('Unstiffened shell')
    c.set_forces(Nsd=N, Msd=M, Tsd=T, Qsd=Q)
    m = c._CylinderMain
    r, t = 6.5, 0.019                                                           # m
    # (2.2.2), (2.2.6), (2.2.7) exactly; (2.2.3) M/(pi r^2 t) at mid-surface, the code uses the outer fibre.
    assert m.sasd == pytest.approx(N * 1e3 / (2 * math.pi * r * t), rel=1e-9)
    assert m.tTsd == pytest.approx(T * 1e3 / (2 * math.pi * r ** 2 * t), rel=1e-9)
    assert m.tQsd == pytest.approx(Q * 1e3 / (math.pi * r * t), rel=1e-9)
    I_tube = math.pi / 4 * ((r + t / 2) ** 4 - (r - t / 2) ** 4)
    assert m.smsd == pytest.approx(M * 1e3 * (r + t / 2) / I_tube, rel=1e-9)
    assert m.smsd == pytest.approx(M * 1e3 / (math.pi * r ** 2 * t), rel=2e-3)


def test_c7_stresses_to_forces_in_kn_and_knm():
    r, t = 6.5, 0.019
    sa, sm, tT, tQ = -50e6, 10e6, 5e6, 2e6                                      # Pa
    N, M, T, Q, _ = helper_cylinder_stress_to_force_to_stress(
        stresses=(sa, sm, tT, tQ, 0), geometry=1, shell_t=t, shell_radius=r, shell_spacing=None,
        CylinderAndCurvedPlate=CylinderAndCurvedPlate)
    assert N == pytest.approx(sa * 2 * math.pi * r * t / 1e3, rel=1e-9)
    assert T == pytest.approx(tT * 2 * math.pi * r ** 2 * t / 1e3, rel=1e-9)
    assert Q == pytest.approx(tQ * math.pi * r * t / 1e3, rel=1e-9)
    assert M == pytest.approx(sm * math.pi * r ** 2 * t / 1e3, rel=2e-3)


# ---------------------------------------------------------------------------------------------------------
# C8 - ring-frame hoop stress index
# ---------------------------------------------------------------------------------------------------------
def test_c8_ring_stress_index_for_ring_stiffened_shell():
    c = build('Ring Stiffened shell', ring=RING)
    c.set_stresses(sasd=-20, tQsd=5, psd=-0.3)
    shRsd = c._CylinderMain.shell_buckling()['shRsd']

    r, t, l = 6500.0, 19.0, 3300.0
    alpha, _, _ = ring_parameters(r, t, l, area_ring=260 * 12 + 49 * 28)
    rf = r - t / 2 - (260 + 28)
    expected = (-0.3 * r / t - NU * -20) / (1 + alpha) * (r / rf)              # (2.2.15) at r_r = r_f
    assert len(shRsd) == 2
    assert shRsd[0] == pytest.approx(expected, rel=1e-9)
    assert math.isnan(shRsd[1])


# ---------------------------------------------------------------------------------------------------------
# C9 - material factor for gammaM0 other than 1.15 / 1.10
# ---------------------------------------------------------------------------------------------------------
def test_c9_material_factor_scales_whole_curve():
    mf = 1.25
    c = build('Unstiffened shell', r=2000, t=20, l=3000, L=3000, mf=mf)
    c.set_stresses(sasd=-100)
    hand_uf, lam, gm = hand_cylinder_shell(2000, 20, 3000, sa=-100, mf=mf)
    assert 0.5 < lam < 1.0
    assert gm == pytest.approx((0.85 + 0.6 * lam) * mf / 1.15, rel=1e-12)
    assert c.get_buckling_results()['Unstiffened shell'] == pytest.approx(hand_uf, rel=1e-9)


# ---------------------------------------------------------------------------------------------------------
# C10 - ring flange eccentricity (3.10.5)
# ---------------------------------------------------------------------------------------------------------
def test_c10_ring_flange_eccentricity_is_checked_dimensionless():
    r, t, l = 1000.0, 10.0, 1000.0
    hw, tw, b, tf = 150.0, 8.0, 100.0, 15.0
    c = build('Ring Stiffened shell', r=r, t=t, l=l, L=10000, ring=(hw, tw, b, tf, 'L', l))
    c.set_stresses(sasd=-20)

    rf = r - t / 2 - (hw + tf)
    Aw, Af = hw * tw, b * tf
    assert hw <= 1.35 * tw * math.sqrt(E / FY)                                 # (3.10.2)
    assert b - tw <= 0.4 * tf * math.sqrt(E / FY)                              # (3.10.3)
    assert hw / tw <= 2 / 3 * math.sqrt(rf * Aw * E / (hw * Af * FY))          # (3.10.4)
    assert (b / 2 - tw / 2) / tw > 1 / 3 * rf / hw * Aw / Af                   # (3.10.5) violated
    assert c.get_buckling_results()['Stiffener check']['ring stiffener'] is False


# ---------------------------------------------------------------------------------------------------------
# C12 - unstiffened curved panel width
# ---------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('s', [680, 1000, 2000])
def test_c12_unstiffened_panel_uses_panel_width(s):
    c = build('Unstiffened panel', panel=s)
    c.set_stresses(sasd=-137.557, tQsd=78.2986, shsd=-3.8934)
    hand_uf = hand_curved_panel_uf(6500, 19, 3300, s, sa=-137.557, tau=78.2986, sh=-3.8934)
    assert c.get_buckling_results()['Unstiffened shell'] == pytest.approx(hand_uf, rel=1e-9)


# ---------------------------------------------------------------------------------------------------------
# C13 - end-cap pressure setting (Table 3-2, notes 1 and 2)
# ---------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('setting, included', [(True, True), (False, False),
                                               ('included in axial stresses', True),
                                               ('not included in axial stresses', False)])
def test_c13_end_cap_setting_selects_pressure_coefficient(setting, included):
    c = build('Unstiffened shell')
    c._CylinderMain.end_cap_pressure_included = setting
    c.set_stresses(psd=-0.15)
    hand_uf, _, _ = hand_cylinder_shell(6500, 19, 3300, sh=-0.15 * 6500 / 19, end_cap_in_axial_stress=included)
    assert c.get_buckling_results()['Unstiffened shell'] == pytest.approx(hand_uf, rel=1e-9)


# ---------------------------------------------------------------------------------------------------------
# C14 - fabrication method of ring frames (3.5.2.6)
# ---------------------------------------------------------------------------------------------------------
def test_c14_fabrication_method_strings_follow_fr_rule():
    def uf(method):
        c = build('Orthogonally Stiffened shell', long=LONG, girder=GIRDER)
        c._CylinderMain.fab_method_ring_stf = method
        c._CylinderMain.fab_method_ring_girder = method
        c.set_stresses(sasd=-60, psd=-0.2)
        return c.get_buckling_results()['Heavy ring frame']

    # Legacy integer code 1 = fabricated (fr = fT); any other code = cold formed (fr = 0.9 fT).
    assert uf('Fabricated') == pytest.approx(uf(1), rel=1e-12)
    assert uf('Cold formed') == pytest.approx(uf(2), rel=1e-12)
    assert uf('Fabricated') < uf('Cold formed')


# ---------------------------------------------------------------------------------------------------------
# C18 - effective shell width of rings (3.5.3)/(3.5.4)
# ---------------------------------------------------------------------------------------------------------
def test_c18_ring_effective_width_limited_to_ring_spacing():
    r, t, l = 6500.0, 19.0, 400.0
    hw, tw, b, tf = 300.0, 12.0, 120.0, 20.0
    c = build('Ring Stiffened shell', r=r, t=t, l=l, ring=(hw, tw, b, tf, 'T', l))
    c.set_stresses(sasd=-40)

    lef = min(1.56 * math.sqrt(r * t) / (1 + 12 * t / r), l)                  # (3.5.3), (3.5.4)
    assert lef == l
    parts = [(lef * t, t / 2, lef * t ** 3 / 12),                              # (area, centroid, own I)
             (hw * tw, t + hw / 2, tw * hw ** 3 / 12),
             (b * tf, t + hw + tf / 2, b * tf ** 3 / 12)]
    area = sum(p[0] for p in parts)
    centroid = sum(p[0] * p[1] for p in parts) / area
    I_ring = sum(p[2] + p[0] * (p[1] - centroid) ** 2 for p in parts)
    assert c._CylinderMain.shell_buckling()['cross section data'][1][4] == pytest.approx(I_ring, rel=1e-9)
