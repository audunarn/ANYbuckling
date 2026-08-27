'''
ANYbuckling - buckling of stiffened flat panels, cylinders and curved
plates according to DNV standards.

Calculation methods:

* Prescriptive flat-plate buckling  - DNV-RP-C201 (:class:`FlatStru`)
* Prescriptive cylinder buckling    - DNV-RP-C202 (:class:`CylStru`)
* Semi-analytical S3/U3 panel solver (:mod:`anybuckling.semianalytical`)
* Optional runner for the DNV PULS Excel sheet (:mod:`anybuckling.puls`)

Quick start::

    from anybuckling import FlatStru

    panel = FlatStru('Flat plate, stiffened')
    panel.set_material(mat_yield=355, emodule=210000, material_factor=1.15)
    panel.set_plate_geometry(spacing=680, thickness=12, span=3300)
    panel.set_stresses(pressure=0.01, sigma_x1=50, sigma_y1=100, tau_xy=5)
    panel.set_stiffener(hw=260, tw=12, bf=49, tf=28, stf_type='bulb', spacing=680)
    panel.set_buckling_parameters(calculation_method='DNV-RP-C201 - prescriptive')
    results = panel.get_buckling_results()
'''
from .api import CylStru, FlatStru
from .helpers import (
    BUCKLING_ACCEPTANCE_TYPES,
    BUCKLING_CALCULATION_METHODS,
    CYLINDER_STRUCTURE_DOMAINS,
    FLAT_STRUCTURE_DOMAINS,
)
from .prescriptive.plates import AllStructure, CalcScantlings, Structure
from .prescriptive.cylinders import CylinderAndCurvedPlate, Shell

__version__ = "0.1.1"

__all__ = [
    "AllStructure",
    "BUCKLING_ACCEPTANCE_TYPES",
    "BUCKLING_CALCULATION_METHODS",
    "CYLINDER_STRUCTURE_DOMAINS",
    "CalcScantlings",
    "CylStru",
    "CylinderAndCurvedPlate",
    "FLAT_STRUCTURE_DOMAINS",
    "FlatStru",
    "Shell",
    "Structure",
]
