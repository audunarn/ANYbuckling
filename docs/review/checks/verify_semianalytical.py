"""Independent checks of anybuckling.semianalytical against classical plate theory.

Run from the repository root after ``pip install -e .``::

    python docs/review/checks/verify_semianalytical.py

Configuration changes are made through ``dataclasses.replace`` on S3SolverConfig;
the package itself is not modified.
"""
import dataclasses
import math

from anybuckling.semianalytical import solver as S

E, NU = 210000.0, 0.3


def header(txt):
    print('\n' + '=' * 100 + '\n' + txt + '\n' + '=' * 100)


def u3_row(a, b, t, sx1, sx2=None, rot=('SS', 'SS')):
    return {'Plate length': a, 'Plate width': b, 'Plate thick.': t, 'Yield stress plate': 355.0,
            'Axial stress': sx1, 'Axial stress 2': sx1 if sx2 is None else sx2,
            'Trans. stress': 0.0, 'Trans. Stress 2': 0.0, 'Shear stress': 0.0,
            'Pressure (fixed)': 0.0, 'In-plane support': 'Integrated',
            'Rotational support': rot[0], 'Rotational support 2': rot[1]}


S3_ROW = {'Length of panel': 3300.0, 'Stiffener spacing': 680.0, 'Plate thick.': 12.0, 'Stiffener type': 'L-bulb',
          'Stiffener boundary': 'Cont', 'Stiff. Height': 260.0, 'Web thick.': 12.0, 'Flange width': 49.0,
          'Flange thick.': 28.0, 'Yield stress plate': 355.0, 'Yield stress stiffener': 355.0,
          'Axial stress': 50.0, 'Trans. stress 1': 100.0, 'Trans. stress 2': 100.0, 'Shear stress': 5.0,
          'Pressure (fixed)': 0.01, 'In-plane support': 'Integrated'}


def sigma_cr_ss(a, b, t):
    """Classical SS plate under uniaxial compression: k = min_m (m b/a + a/(m b))^2."""
    D = E * t ** 3 / (12 * (1 - NU ** 2))
    k = min((m * b / a + a / (m * b)) ** 2 for m in range(1, 80))
    return k * math.pi ** 2 * D / (b ** 2 * t)


header('S1  Ritz basis fixed at m<=3, n<=2: U3 elastic factor vs closed form (b=1000, t=12, sigma_x=100)')
rich = dataclasses.replace(S.S3SolverConfig(), longitudinal_modes=tuple(range(1, 13)))
for ab in (1, 3, 4, 5, 8, 12):
    a = 1000.0 * ab
    panel = S.row_to_u3_input(u3_row(a, 1000.0, 12.0, 100.0))
    vals = []
    for cfg in (S.S3SolverConfig(), rich):
        f = S.elastic_u3_buckling_factors(panel, S.build_u3_ritz_modes(panel, cfg), cfg)
        vals.append(f['critical_factor'])
    ref = sigma_cr_ss(a, 1000.0, 12.0) / 100.0
    print(f'  a/b={ab:2d}: closed form={ref:.4f}  default basis={vals[0]:.4f} ({100 * (vals[0] / ref - 1):+6.1f}%)  '
          f'm<=12 basis={vals[1]:.4f} ({100 * (vals[1] / ref - 1):+.1f}%)')

header('S2  Stress gradients replaced by their mean: in-plane bending psi=-1 gives no buckling driver')
panel = S.row_to_u3_input(u3_row(1000.0, 1000.0, 12.0, 100.0, -100.0))
cfg = S.S3SolverConfig()
f = S.elastic_u3_buckling_factors(panel, S.build_u3_ritz_modes(panel, cfg), cfg)
D = E * 12.0 ** 3 / (12 * (1 - NU ** 2))
print(f'  sigma_x1=+100, sigma_x2=-100: stress seen by solver={panel.axial_stress}, critical factor={f["critical_factor"]}; '
      f'classical k~23.9 -> factor {23.9 * math.pi ** 2 * D / (1000.0 ** 2 * 12.0) / 100:.2f}')

header('S3  U3 rotational supports accepted but not used (SS vs CL identical)')
for rot in (('SS', 'SS'), ('CL', 'CL')):
    r = S.solve_u3_panel(S.row_to_u3_input(u3_row(3000.0, 1000.0, 12.0, 100.0, rot=rot)))
    print(f'  {rot}: buckling UF={r.buckling_usage_factor:.6f} ultimate UF={r.ultimate_usage_factor:.6f}')

header('S4  README S3 panel: sensitivity to mode set and the 1.10 local-plate restraint constant')
cases = [('default (1,2,3)x(1,2), restraint 1.10', S.S3SolverConfig()),
         ('(1..7)x(1,2,3), restraint 1.10', dataclasses.replace(S.S3SolverConfig(), longitudinal_modes=tuple(range(1, 8)),
                                                               transverse_modes=(1, 2, 3))),
         ('(1..7)x(1,2,3), restraint 1.00', dataclasses.replace(S.S3SolverConfig(), longitudinal_modes=tuple(range(1, 8)),
                                                               transverse_modes=(1, 2, 3),
                                                               flanged_local_plate_restraint_factor=1.0))]
for sx, sy in ((50.0, 100.0), (150.0, 30.0)):
    row = dict(S3_ROW, **{'Axial stress': sx, 'Trans. stress 1': sy, 'Trans. stress 2': sy})
    for label, c in cases:
        r = S.solve_s3_panel(S.row_to_s3_input(row), c)
        print(f'  sx={sx:5.1f} sy={sy:5.1f}  {label:40s} buckling UF={r.buckling_usage_factor:.4f} '
              f'ultimate UF={r.ultimate_usage_factor:.4f}')

header('S5  effective_stiffener_plate_width uses r**1.6; IACS/CSR form is 1.12/(1+1.75 r**-1.6)')
for length in (2355.0, 3400.0, 5889.0):
    panel = S.row_to_s3_input(dict(S3_ROW, **{'Length of panel': length}))
    out = S.effective_stiffener_plate_width(panel)
    r = out['length_ratio']
    ref = min(1.0, 1.12 / (1 + 1.75 * r ** -1.6))
    print(f'  L={length:6.0f} (l_eff/s={r:4.2f}): code coefficient={out["coefficient"]:.3f}  CSR={ref:.3f}  '
          f'(only active with use_effective_stiffener_width=True)')


# ------------------------------------------------------------------------------------------------
# Ultimate path, adapters and API (through FlatStru, as a user would call it)
# ------------------------------------------------------------------------------------------------
from anybuckling import FlatStru


def flat(s, t, l, stf, sx=0.0, sx2=None, sy=0.0, tau=0.0, p=0.0, method='SemiAnalytical S3/U3',
         stiffener_support='Continuous', **bp):
    o = FlatStru('Flat plate, stiffened')
    o.set_material(355, 210000, 1.15, 0.3)
    o.set_plate_geometry(spacing=s, thickness=t, span=l)
    o.set_stresses(pressure=p, sigma_x1=sx, sigma_x2=sx if sx2 is None else sx2, sigma_y1=sy, sigma_y2=sy,
                   tau_xy=tau)
    hw, tw, bf, tf, ty = stf
    o.set_stiffener(hw=hw, tw=tw, bf=bf, tf=tf, stf_type=ty, spacing=s)
    o.set_buckling_parameters(calculation_method=method, buckling_acceptance='ultimate',
                              stiffener_support=stiffener_support, **bp)
    if method.startswith('Semi'):
        o.set_puls_parameters(stiffener_end=stiffener_support)
    return o


def euler_full_plate(s, t, L, stf):
    hw, tw, bf, tf, _ = stf
    parts = [(s * t, 0.0, s * t ** 3 / 12), (hw * tw, t / 2 + hw / 2, tw * hw ** 3 / 12),
             (bf * tf, t / 2 + hw + tf / 2, bf * tf ** 3 / 12)]
    a = sum(p[0] for p in parts)
    zc = sum(p[0] * p[1] for p in parts) / a
    inertia = sum(i + ar * (z - zc) ** 2 for ar, z, i in parts)
    return math.pi ** 2 * E * inertia / (a * L ** 2)


header('U1  Ultimate capacity above the column Euler load (solver.py:4913-4925, 5015-5020)')
for stf, label in (((120, 10, 0, 0, 'FB'), 'FB 120x10'), ((120, 14, 80, 12, 'T'), 'T 120x14/80x12')):
    r = flat(800, 16, 7000, stf, sx=100).get_buckling_results()
    sE = euler_full_plate(800.0, 16.0, 7000.0, stf)
    c201 = flat(800, 16, 7000, stf, sx=100, method='DNV-RP-C201 - prescriptive').get_buckling_results()['Stiffener']
    print(f'  {label:15s} L=7 m, sigma_x=100: selected (ultimate) UF={r["selected UF"]:.3f} '
          f'(raw {r["ultimate UF raw"]:.4f}, valid={r["valid prediction"]}, confidence={r["confidence"]}); '
          f'Euler stress with full plate={sE:.1f} MPa -> UF >= {100 / sE:.2f}; RP-C201 UF={c201["Overpressure plate side"]:.2f}')

header('U2  Pressure-dominated rule lifts the ultimate capacity to the elastic factor (solver.py:5809-5839)')
for p in (0.0, 0.05, 0.10, 0.40):
    r = flat(500, 18, 1500, (200, 12, 100, 15, 'T'), sx=200, p=p).get_buckling_results()
    print(f'  T 200x12/100x15, sigma_x=200, p={p:4.2f} MPa: ultimate UF raw={r["ultimate UF raw"]:.4f} '
          f'(sigma_x/fy alone = {200 / 355:.3f})')

header('U3  Load factors are ignored by the SemiAnalytical path (api.py:164-217, solver.py:927-998)')
for lf in (1.0, 1.5):
    kw = dict(s=680, t=12, l=3300, stf=(260, 12, 49, 28, 'bulb'), sx=50, sy=100, tau=5, p=0.01,
              load_factor_stresses=lf, load_factor_pressure=lf)
    r = flat(**kw).get_buckling_results()
    c = flat(method='DNV-RP-C201 - prescriptive', **kw).get_buckling_results()
    print(f'  load factors={lf}: SemiAnalytical ultimate UF={r["ultimate UF"]:.4f} | RP-C201 plate UF={c["Plate"]["Plate buckling"]:.4f}')

header('U4  Sniped stiffener: ultimate UF falls as lateral pressure rises (plate-induced branch sign, solver.py:5008-5014)')
for p in (0.0, 0.05, 0.10):
    r = flat(800, 10, 3500, (200, 12, 250, 20, 'T'), sx=100, p=p, stiffener_support='Sniped').get_buckling_results()
    print(f'  T 200x12/250x20 sniped, sigma_x=100, p={p:4.2f}: ultimate UF raw={r["ultimate UF raw"]:.4f}')

header('U5  Mixed-sign axial stress keeps only the compression end (solver.py:921-924)')
for sx1, sx2 in ((20, -300), (20, 20)):
    r = flat(680, 12, 3300, (260, 12, 49, 28, 'bulb'), sx=sx1, sx2=sx2).get_buckling_results()
    print(f'  sigma_x1={sx1}, sigma_x2={sx2}: ultimate UF raw={r["ultimate UF raw"]:.4f} (300 MPa tension alone >= {300 / 355:.3f})')

header('U6  CSR proportion check reads thicknesses < 10 mm as metres (solver.py:890-896)')
r = flat(700, 8, 3300, (400, 8, 120, 8, 'T'), sx=50, p=0.01).get_buckling_results()
v = r['CSR requirement']['values']
print(f'  t=tw=tf=8 mm: plate net thickness used={v["plate_net_thickness"]:.0f} mm, CSR vector={r["CSR"]} ({r["CSR color"]}); '
      f'hand: s/t_net = 700/6 = {700 / 6:.1f} > limit {v["plate_slenderness_limit"]:.1f} -> should fail')
