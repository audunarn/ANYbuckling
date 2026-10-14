"""Independent verification of anybuckling.prescriptive.cylinders (DNV-RP-C202) via the public CylStru API.

Run from the repository root after ``pip install -e .``::

    python docs/review/checks/verify_c202_cylinders.py

Each block contrasts two inputs that the standard says must give a different (or the
same) answer, or compares with a hand calculation. All blocks show the behaviour after the
fixes on this branch (the review report lists the values before the fixes). Nothing in the
package is modified.
"""
import math

from anybuckling import CylStru
from anybuckling.prescriptive.cylinders import CylinderAndCurvedPlate


def cyl(domain='Unstiffened shell', r=6500, t=19, l=3300, L=20000, long=None, ring=None, girder=None,
        panel=680, mf=1.15, fab=('Fabricated', 'Fabricated'), endcap=True):
    c = CylStru(calculation_domain=domain)
    c.set_material(mat_yield=355, emodule=210000, material_factor=mf, poisson=0.3)
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
    c.set_shell_buckling_parmeters()
    return c


def uf(c, key):
    return c.get_buckling_results()[key]


def header(txt):
    print('\n' + '=' * 100 + '\n' + txt + '\n' + '=' * 100)


LONG = (232, 12, 49, 28, 'bulb', 680.367)
RING = (260, 12, 49, 28, 'bulb', 680)
GIRDER = (500, 15, 200, 25, 'T', 3300)

header('C1  (fixed) Ring-stiffened shell includes the additional hoop stress shsd, like the unstiffened shell')
for shsd in (0, -50):
    u = cyl('Unstiffened shell')
    u.set_stresses(sasd=-40, tQsd=10, shsd=shsd)
    rs = cyl('Ring Stiffened shell', ring=RING)
    rs.set_stresses(sasd=-40, tQsd=10, shsd=shsd)
    print(f'  shsd={shsd:4d} MPa: unstiffened shell UF={uf(u, "Unstiffened shell"):.3f} | '
          f'ring-stiffened, shell-between-rings UF={uf(rs, "Unstiffened shell"):.3f}')

header('C2  (fixed) Conical shell, stress input: the additional hoop stress shsd is kept')
for shsd in (0, -50, -100):
    c = CylStru('Unstiffened conical shell')
    c.set_material(355, 210000, 1.15, 0.3)
    c.set_imperfection()
    c.set_fabrication_method()
    c.set_end_cap_pressure_included_in_stress()
    c.set_uls_or_als()
    c.set_conical_shell_geometry(r1=4000, r2=6500, length=5000, thickness=20)
    c.set_conical_stresses(sasd=-60, shsd=shsd)
    print(f'  shsd={shsd:4d} MPa: cone UF={uf(c, "Unstiffened conical shell"):.4f}')

header('C3  (fixed) Longitudinally stiffened shell is continuous at sigma_x = 0')
for sa in (0, -0.5, -5):
    c = cyl('Longitudinal Stiffened shell', long=LONG)
    c.set_stresses(sasd=sa, tQsd=50, shsd=-100)
    print(f'  sasd={sa:5}: UF={uf(c, "Longitudinal stiffened shell"):.3f}')

header('C4  (fixed) Longitudinally stiffened full shell is independent of the sign of the bending stress')
for sm in (50, -50):
    c = cyl('Longitudinal Stiffened shell', long=LONG)
    c.set_stresses(sasd=-20, smsd=sm, shsd=-100)
    print(f'  smsd={sm:4d}: UF={uf(c, "Longitudinal stiffened shell"):.3f}  (a full cylinder sees +/-sm around the circumference)')

header('C5  (fixed) Unstiffened shells get the Sec 3.8 column buckling check')
c = cyl('Unstiffened shell', r=500, t=12, l=3000, L=40000)
c.set_stresses(sasd=-150)
res = c.get_buckling_results()
ic = 500 / math.sqrt(2)
print(f'  tube r=500 t=12 L=40 m, sigma_a=150: shell UF={res["Unstiffened shell"]:.3f}, '
      f'column UF={res["Column stability UF"]}, column check={res["Column stability check"]}')
print(f'  hand: (L/ic)^2={(40000 / ic) ** 2:.0f} >= 2.5E/fy={2.5 * 210000 / 355:.0f} -> check required; '
      f'Euler stress={math.pi ** 2 * 210000 / (40000 / ic) ** 2:.1f} MPa')

header('C6  (fixed) Column curve for lambda > 1.34: 0.9 fak/lambda^2, eq. (3.8.6)')
for lam in (1.34, 1.3401, 1.5, 2.0):
    code = CylinderAndCurvedPlate._column_reduced_buckling_strength(1.0, lam)
    ref = (1 - 0.28 * lam ** 2) if lam <= 1.34 else 0.9 / lam ** 2
    print(f'  lambda={lam:<6}: code fkc/fak={code:.4f}  RP={ref:.4f}  ({100 * (code / ref - 1):+.1f}%)')

header('C7  (fixed) set_forces: bending and torsion stresses in the right units')
c = cyl('Unstiffened shell')
c.set_forces(Nsd=0, Msd=10000, Tsd=10000, Qsd=0, psd=0)
m = c._CylinderMain
r_, t_ = 6500.0, 19.0
print(f'  M=T=10000 kNm, r=6.5 m, t=19 mm: code sigma_m={m.smsd / 1e6:.3f} MPa (thin-wall hand {10000e6 * r_ / (math.pi * r_ ** 3 * t_):.3f}; code uses r + t/2); '
      f'code tau_T={m.tTsd / 1e6:.3f} MPa (hand {10000e6 / (2 * math.pi * r_ ** 2 * t_):.3f})')

header('C12 (fixed) Unstiffened curved panel uses the panel width s (reference: 0.855, 1.292, 2.625)')
for sp in (680, 1000, 2000):
    c = cyl('Unstiffened panel', panel=sp)
    c.set_stresses(sasd=-137.557, tQsd=78.2986, shsd=-3.8934)
    print(f'  panel spacing={sp:5d}: UF={uf(c, "Unstiffened shell"):.4f}')

header('C13/C14 (fixed) End-cap flag (True: psi = 4, False: psi = 2) and fabrication method (fT vs 0.9 fT) take effect')
for ec in (True, False):
    c = cyl('Unstiffened shell', endcap=ec)
    c.set_stresses(psd=-0.15)
    print(f'  end cap pressure included={ec!s:5}: UF={uf(c, "Unstiffened shell"):.4f}')
for fab in (('Fabricated', 'Fabricated'), ('Cold formed', 'Cold formed')):
    c = cyl('Orthogonally Stiffened shell', long=LONG, girder=GIRDER, fab=fab)
    c.set_stresses(sasd=-60, psd=-0.2)
    print(f'  fabrication={fab[0]:12s}: heavy ring frame UF={uf(c, "Heavy ring frame"):.4f}')
