"""Independent verification of anybuckling.prescriptive.plates against DNV-RP-C201 (Oct 2010).

Run from the repository root after ``pip install -e .``::

    python docs/review/checks/verify_c201_plates.py

Each block prints the code value next to an independent hand calculation or a
targeted counter-example. Nothing in the package is modified; where a corrected
formula is evaluated, it is evaluated here, not patched into the package
(except D3, which temporarily swaps one method and restores it).
"""
import inspect
import math

from scipy.optimize import minimize
import scipy.optimize as so

from anybuckling import FlatStru
from anybuckling.prescriptive import plates as P

E, NU, FY, GM = 210000.0, 0.3, 355.0, 1.15
G_MOD = E / (2 * (1 + NU))


def build(domain='Flat plate, stiffened', s=680, t=12, l=3300, sx=50, sy=100, tau=5, p=0.01,
          stf=(260, 12, 49, 28, 'bulb'), gir=(500, 15, 200, 25, 'T'), **bp):
    o = FlatStru(domain)
    o.set_material(mat_yield=FY, emodule=E, material_factor=GM, poisson=NU)
    o.set_plate_geometry(spacing=s, thickness=t, span=l)
    o.set_stresses(pressure=p, sigma_x1=sx, sigma_x2=sx, sigma_y1=sy, sigma_y2=sy, tau_xy=tau)
    if domain != 'Flat plate, unstiffened':
        hw, tw, bf, tf, ty = stf
        o.set_stiffener(hw=hw, tw=tw, bf=bf, tf=tf, stf_type=ty, spacing=s)
    if domain == 'Flat plate, stiffened with girder':
        hw, tw, bf, tf, ty = gir
        o.set_girder(hw=hw, tw=tw, bf=bf, tf=tf, stf_type=ty, spacing=700)
    kw = dict(calculation_method='DNV-RP-C201 - prescriptive', buckling_acceptance='ultimate')
    kw.update(bp)
    o.set_buckling_parameters(**kw)
    return o


def instrumented(method_name):
    """Copy of an AllStructure method that also returns its locals (for inspection only)."""
    lines = inspect.getsource(getattr(P.AllStructure, method_name)).splitlines()
    src = '\n'.join(ln[4:] if ln.startswith('    ') else ln for ln in lines) + '\n'
    ret = '    return girder_data\n' if method_name == 'girder' else '    return stf_pl_data\n'
    src = src.replace(ret, ret.rstrip('\n') + ', dict(locals())\n')
    ns = dict(P.__dict__)
    ns['minimize'] = minimize
    exec(src.replace(f'def {method_name}(', f'def {method_name}_dbg('), ns)
    return ns[f'{method_name}_dbg']


def header(txt):
    print('\n' + '=' * 100 + '\n' + txt + '\n' + '=' * 100)


# --------------------------------------------------------------------------------------------
header('V0  Hand calculation of the README stiffened panel (bulb 260x12/49x28) - intermediate values')
o = build()
A = o._FlatStructure
up = A.unstiffened_plate_buckling()
sp = A.stiffened_panel(unstf_pl_data=up)
s, t, l = 680.0, 12.0, 3300.0
hw, tw, bf, tf = 260.0, 12.0, 49.0, 28.0
sx, sy, tau = 50.0, 100.0, 5.0
lp = 0.525 * s / t * math.sqrt(FY / E)
Cxs = (lp - 0.22) / lp ** 2 if lp > 0.673 else 1.0
lc = 1.1 * s / t * math.sqrt(FY / E)
kap = 1 / (2 * lc ** 2) + 0.07 if lc >= 2 else None
syR = (1.3 * t / l * math.sqrt(E / FY) + kap * (1 - 1.3 * t / l * math.sqrt(E / FY))) * FY
ci = 1 - s / (120 * t)
arg = 1 - (sy / syR) ** 2 + ci * sx * sy / (Cxs * FY * syR)
Cys = math.sqrt(max(arg, 0.0))
def section(be):
    parts = [(be * t, 0.0, be * t ** 3 / 12), (hw * tw, t / 2 + hw / 2, tw * hw ** 3 / 12),
             (bf * tf, t / 2 + hw + tf / 2, bf * tf ** 3 / 12)]
    a = sum(p[0] for p in parts)
    zc = sum(p[0] * p[1] for p in parts) / a
    return a, zc, sum(i + ar * (z - zc) ** 2 for ar, z, i in parts)
Is = section(s)[2]
ef = bf / 2 - tw / 2
Ipo = tw * hw ** 3 / 3 + hw * tw ** 3 / 12 + bf * tf * (hw + tf / 2) ** 2 + bf * tf ** 3 / 12 \
    + tf * ((ef + bf / 2) ** 3 - (ef - bf / 2) ** 3) / 3
Aw, Af = hw * tw, bf * tf
yc = (Aw * tw / 2 + Af * bf / 2) / (Aw + Af)
Iz = hw * tw ** 3 / 12 + Aw * (tw / 2 - yc) ** 2 + tf * bf ** 3 / 12 + Af * (bf / 2 - yc) ** 2
for name, ref, code in [('lambda_p (7.3)', lp, sp['alphap']), ('Cxs (7.12)', Cxs, sp['Cxs']),
                        ('sigma_yR (6.8)', syR, up['syR']), ('Cys (7.14, clamped >= 0)', Cys, sp['Cys_comp']),
                        ('Is full plate [mm4]', Is, sp['Is']),
                        ('Ipo [mm4]', Ipo, A._Stiffener.get_polar_moment()),
                        ('Iz [mm4]', Iz, A._Stiffener.get_Iz_moment_of_inertia())]:
    ok = abs(ref - code) <= 1e-6 * max(1.0, abs(ref))
    print(f'  {name:28s} hand={ref:14.6g}  code={code:14.6g}  {"MATCH" if ok else "DIFF"}')
print(f'  note: Cys argument = {arg:.3f} < 0, i.e. sigma_y = 100 MPa exceeds sigma_yR = {syR:.1f} MPa; '
      f'effective plate width is zero for this README example.')

# --------------------------------------------------------------------------------------------
header('D1  Sec 6.5: C_tau_e is computed but overwritten (plates.py:1861-1862); combined check always uses C_tau')
o = build('Flat plate, unstiffened', s=1000, t=10, l=3000, sx=30, sy=20, tau=40, p=0)
d = o._FlatStructure.unstiffened_plate_buckling()
s, t, l = 1000.0, 10.0, 3000.0
kl = 5.34 + 4 * (s / l) ** 2
lw = 0.795 * s / t * math.sqrt(FY / E / kl)
Cte = 1 if lw <= 0.8 else (1 - 0.8 * (lw - 0.8) if lw <= 1.25 else 1 / lw ** 2)
sxRd = 30 / d['UF Longitudinal stress']
syRd = d['syRd']
ci = max(0.0, 1 - s / (120 * t))
uf_e = (30 / sxRd) ** 2 + (20 / syRd) ** 2 - ci * (30 / sxRd) * (20 / syRd) \
    + (40 / (Cte * FY / GM / math.sqrt(3))) ** 2
print(f'  s/t=100, sigma_y compressive, lambda_w={lw:.3f}: UF combined code={d["UF Combined stresses"]:.3f}, '
      f'with C_tau_e={uf_e:.3f}  ({100 * (uf_e / d["UF Combined stresses"] - 1):+.0f}%)')

# --------------------------------------------------------------------------------------------
header('D2  Sec 6.4: C_tau middle branch coded as 1-0.675(lw-0.8) up to 1.25; continuous form is 1-0.625(lw-0.8) up to 1.2')
for lw in (1.0, 1.2, 1.22, 1.25):
    code = 1 - 0.675 * (lw - 0.8) if lw < 1.25 else 0.9 / lw
    ref = 1 - 0.625 * (lw - 0.8) if lw <= 1.2 else 0.9 / lw
    print(f'  lambda_w={lw:.2f}: code={code:.4f}  ref={ref:.4f}  ({100 * (code / ref - 1):+.1f}%, conservative)')
print('  continuity check of code form at 1.25: left=%.4f right=%.4f' % (1 - 0.675 * 0.45, 0.9 / 1.25))

# --------------------------------------------------------------------------------------------
header('D3  Sec 7.5.2: St. Venant It omits (1 - 0.63 tf/bf) on the flange term (plates.py:651)')
orig = P.Structure.get_torsional_moment_venant
def It_rp(self, reduced_tw=None, efficient_flange=True):
    tf_ = self._flange_th * 1000
    tw_ = self._web_th * 1000 if reduced_tw is None else reduced_tw
    bf_ = self._flange_width * 1000
    hw_ = self._web_height * 1000
    return (hw_ * tw_ ** 3 + (0 if bf_ == 0 else bf_ * tf_ ** 3 * (1 - 0.63 * tf_ / bf_))) / 3
try:
    for label, kw in [('README bulb, sigma_y=100', {}), ('README bulb, sigma_y=60', dict(sy=60)),
                      ('T 400x12/150x20', dict(s=750, t=15, sx=150, sy=30, p=0.05, stf=(400, 12, 150, 20, 'T')))]:
        P.Structure.get_torsional_moment_venant = orig
        r0 = build(**kw).get_buckling_results()['Stiffener']
        P.Structure.get_torsional_moment_venant = It_rp
        r1 = build(**kw).get_buckling_results()['Stiffener']
        print(f'  {label:28s} stiffener-side UF code={r0["Overpressure stiffener side"]:.3f} -> with RP It '
              f'{r1["Overpressure stiffener side"]:.3f}; plate side {r0["Overpressure plate side"]:.3f} -> '
              f'{r1["Overpressure plate side"]:.3f}')
finally:
    P.Structure.get_torsional_moment_venant = orig

# --------------------------------------------------------------------------------------------
header('D4  Sec 7.7.1: z* optimised separately for the plate-side and stiffener-side sets (plates.py:2201, 2220)')
captured = []
real_minimize = so.minimize
def spy(fun, x0, bounds=None, **k):
    captured.append((fun, bounds))
    return real_minimize(fun, x0, bounds=bounds, **k)
P_minimize_name = 'minimize'
try:
    so.minimize = spy
    for label, kw in [('README', {}),
                      ('T 400x12/150x20', dict(s=750, t=15, sx=150, sy=30, p=0.05, stf=(400, 12, 150, 20, 'T'))),
                      ('FB 200x14', dict(s=600, t=12, sx=120, sy=0, tau=10, p=0.1, stf=(200, 14, 0, 0, 'FB')))]:
        captured.clear()
        r = build(**kw).get_buckling_results()['Stiffener']
        (fpl, b), (fst, _) = captured[-2], captured[-1]
        lo, hi = b[0]
        zs = [lo + (hi - lo) * i / 20000 for i in range(20001)]
        joint = min(max(fpl(z), fst(z)) for z in zs)
        print(f'  {label:16s} code: plate side={r["Overpressure plate side"]:.3f} stiffener side='
              f'{r["Overpressure stiffener side"]:.3f}  | one z* for eqs 7.50-7.57: {joint:.3f}')
finally:
    so.minimize = real_minimize

# --------------------------------------------------------------------------------------------
GD = 'Flat plate, stiffened with girder'
header('D5  API: girder span LG fixed at 10 m (api.py:264,336,376); set_girder(spacing=) has no effect')
for spc in (700, 9000):
    o = build(GD, sy=40)
    o.set_girder(hw=500, tw=15, bf=200, tf=25, stf_type='T', spacing=spc)
    print(f'  set_girder(spacing={spc}) -> {o.get_buckling_results()["Girder"]}')
for lg in (4.0, 10.0, 16.0):
    o = build(GD, sy=40)
    o._FlatStructure.Plate.girder_lg = lg
    print(f'  Plate.girder_lg={lg:4.1f} m (private attribute) -> {o.get_buckling_results()["Girder"]}')

girder_dbg = instrumented('girder')
def girder_locals(**kw):
    lg = kw.pop('lg', None)
    A = build(GD, **kw)._FlatStructure
    if lg is not None:
        A.Plate.girder_lg = lg
    up = A.unstiffened_plate_buckling()
    sp = A.stiffened_panel(unstf_pl_data=up)
    return girder_dbg(A, unstf_pl_data=up, stf_pl_data=sp)

header('D6  Sec 7.7.2 applied to sniped girders: equation-set selection inverted (plates.py:2481-2486)')
for sx_, sy_ in ((0, 40), (5, 30)):
    gd, L = girder_locals(sx=sx_, sy=sy_, tau=0, p=0.0, girder_support='Sniped')
    q8 = L['qSd_girder_side'] * L['Lg'] ** 2 / 8
    nz = L['NySd'] * L['zstar']
    ref = max(0, L['uf_7_60'], L['uf_7_61']) if q8 >= nz else max(0, L['uf_7_62'], L['uf_7_63'])
    print(f'  sx={sx_} sy={sy_}: q LG^2/8={q8:.3e} < NySd z*={nz:.3e}; code UF={gd["UF Simplified girder side"]:.3f}, '
          f'RP selection (7.62/7.63) UF={ref:.3f}')

header('D7  Sec 8.4: shear-lag limit le <= 0.3 LG / 0.4 LG computed (eff_width_sec_mod) but never used')
gd, L = girder_locals(sx=10, sy=10, tau=5, p=0.002, lg=6.0)
Gd = L['self']._Girder
le_lim = L['tot_min_lim']
Iy = Gd.get_moment_of_intertia(efficent_se=le_lim / 1000) * 1e12
zp = Gd.get_cross_section_centroid_with_effective_plate(le_lim / 1000) * 1000 - 6
print(f'  LG=6 m: le used={L["le"]:.0f} mm, 0.3 LG={0.3 * L["Lg"]:.0f} mm; Wep code={L["Wep"]:.3e} mm3 vs '
      f'with limit={Iy / zp:.3e} mm3 ({100 * (L["Wep"] / (Iy / zp) - 1):+.0f}%)')

header('D8  Sec 8.2: p0 for sigma_x in tension uses (l - s/LG) and gammaM instead of (1 - s/LG) and fy (plates.py:2397)')
t, s, l, Lg, hwG, As = 12.0, 680.0, 3300.0, 10000.0, 500.0, 260 * 12 + 49 * 28
ref = 0.4 * (t + As / s) / (hwG * (1 - s / Lg)) * FY / E * (Lg / l) ** 2
code = 0.4 * (t + As / s) / (hwG * (l - s / Lg)) * GM / E * (Lg / l) ** 2
print(f'  coefficient on C*tau_Sd: code={code:.3e}, compression-form={ref:.3e} (code is {ref / code:.1e}x smaller)')

header('D9  Sec 8.2: tau_crg defaults to 0.6 fy when panel length Lp is not given (Lp=0 branch, plates.py:2301-2304)')
for Lp in (None, 6 * 3300, 12 * 3300):
    o = build(GD, sx=60, sy=40, tau=80, panel_length_Lp=Lp)
    print(f'  Lp={Lp}: girder plate-side UF={o.get_buckling_results()["Girder"]["Overpressure plate side"]:.3f}')

header('D10 Sec 7.7: reduced web thickness tw(1-V/VRd) switches on at V/VRd=0.5 (jump) and crashes for V/VRd>=1')
for p in (0.094, 0.0955, 0.19, 0.20):
    try:
        r = build(s=680, t=12, l=3300, sx=50, sy=20, tau=5, p=p, stf=(200, 6, 0, 0, 'FB')).get_buckling_results()
        print(f'  p={p}: V/VRd={r["Stiffener"]["Shear capacity"]:.3f} plate-side UF={r["Stiffener"]["Overpressure plate side"]:.2f}')
    except ZeroDivisionError as err:
        print(f'  p={p}: get_buckling_results() raised ZeroDivisionError ({err})')

header('D11 Special provisions: sigma_y sampled at the low-stress end; Z floor = 15 mm3 instead of 15 cm3')
import copy, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'tests'))
import example_data as ex
dd = copy.deepcopy(ex.obj_dict)
dd['panel or shell'] = ['panel', '']
sc = P.CalcScantlings(dd)
for y1, y2 in ((200, 0), (0, 200)):
    sc.set_stresses(sigy1=y1, sigy2=y2, sigx1=50, sigx2=50, tauxy=0)
    print(f'  sigma_y1={y1:3d}, sigma_y2={y2:3d}: Z_min={sc.get_dnv_min_section_modulus(200) * 1e9:8.0f} mm3, '
          f't_min={sc.get_dnv_min_thickness(200):.2f} mm  (same stress field, different answer)')
print(f'  floor used: {15 ** 3} mm3; DNV-OS-C101 floor: 15000 mm3')
