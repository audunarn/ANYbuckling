# ANYbuckling

Buckling of stiffened flat panels, cylinders and curved plates according
to DNV standards.

ANYbuckling is the standalone calculation engine extracted from
[ANYstructure](https://github.com/audunarn/ANYstructure). It has no GUI
and only depends on numpy and scipy.

## Calculation methods

| Method | Standard | Entry point |
|---|---|---|
| Prescriptive flat-plate buckling | DNV-RP-C201 | `anybuckling.FlatStru` |
| Prescriptive cylinder / curved plate buckling | DNV-RP-C202 | `anybuckling.CylStru` |
| Semi-analytical stiffened/unstiffened panel (S3/U3) | PULS-type ultimate capacity | `anybuckling.semianalytical` |
| DNV PULS Excel runner (optional, licensed sheet required) | PULS | `anybuckling.puls.PULSpanel` |

The ML-based buckling predictions available in ANYstructure are *not*
part of this package.

## Installation

```bash
pip install ANYbuckling
```

Optional extras:

```bash
pip install ANYbuckling[excel]   # xlwings, needed only for the PULS Excel runner
```

## Quick start

### Flat stiffened panel (DNV-RP-C201, prescriptive)

```python
from anybuckling import FlatStru

panel = FlatStru('Flat plate, stiffened')
panel.set_material(mat_yield=355, emodule=210000, material_factor=1.15, poisson=0.3)
panel.set_plate_geometry(spacing=680, thickness=12, span=3300)
panel.set_stresses(pressure=0.01, sigma_x1=50, sigma_x2=50,
                   sigma_y1=100, sigma_y2=100, tau_xy=5)
panel.set_stiffener(hw=260, tw=12, bf=49, tf=28, stf_type='bulb', spacing=680)
panel.set_buckling_parameters(calculation_method='DNV-RP-C201 - prescriptive',
                              buckling_acceptance='ultimate')

results = panel.get_buckling_results()
print(results['Plate']['Plate buckling'])          # plate UF
print(results['Stiffener'])                        # stiffener UFs
```

Domains: `'Flat plate, unstiffened'`, `'Flat plate, stiffened'` and
`'Flat plate, stiffened with girder'`. Switch `calculation_method` to
`'SemiAnalytical S3/U3'` to run the same panel through the
semi-analytical solver.

### Cylinder (DNV-RP-C202)

```python
from anybuckling import CylStru

cyl = CylStru(calculation_domain='Longitudinal Stiffened shell')
cyl.set_material(mat_yield=355, emodule=210000, material_factor=1.15, poisson=0.3)
cyl.set_stresses(sasd=-137.6, tQsd=78.3, shsd=-3.9)
cyl.set_shell_geometry(radius=6500, thickness=19,
                       tot_length_of_shell=20000, distance_between_rings=3300)
cyl.set_longitudinal_stiffener(hw=232, tw=12, bf=49, tf=28, spacing=680)
cyl.set_imperfection()
cyl.set_fabrication_method()
cyl.set_end_cap_pressure_included_in_stress()
cyl.set_uls_or_als('ULS')
cyl.set_shell_buckling_parmeters()

print(cyl.get_buckling_results())
```

Domains cover unstiffened, longitudinally stiffened, ring stiffened and
orthogonally stiffened shells and panels, plus unstiffened conical
shells, with either stress or force input.

### Semi-analytical S3/U3 solver directly

```python
from anybuckling.semianalytical import row_to_s3_input, solve_s3_panel

row = {
    'Length of panel': 3300.0, 'Stiffener spacing': 680.0, 'Plate thick.': 12.0,
    'Stiffener type': 'L-bulb', 'Stiffener boundary': 'Cont',
    'Stiff. Height': 260.0, 'Web thick.': 12.0, 'Flange width': 49.0,
    'Flange thick.': 28.0, 'Yield stress plate': 355.0, 'Yield stress stiffener': 355.0,
    'Axial stress': 50.0, 'Trans. stress 1': 100.0, 'Trans. stress 2': 100.0,
    'Shear stress': 5.0, 'Pressure (fixed)': 0.01, 'In-plane support': 'Integrated',
}
result = solve_s3_panel(row_to_s3_input(row))
print(result.valid, result.buckling_usage_factor, result.ultimate_usage_factor)
```

The solver also ships a benchmark and verification CLI:

```bash
python -m anybuckling.semianalytical --help
```

## Relation to ANYstructure

The classes keep their original ANYstructure names (`Structure`,
`CalcScantlings`, `AllStructure`, `Shell`, `CylinderAndCurvedPlate`,
`PULSpanel`) because ANYstructure uses this package as a dependency. Fatigue,
load combinations, optimization, FE integration, machine-learning prediction,
and the GUI remain in ANYstructure.

## Development

```bash
pip install -e .[dev]
pytest
```

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
