'''Longitudinally stiffened cylindrical shell checked to DNV-RP-C202.'''
from anybuckling import CylStru


if __name__ == '__main__':
    cyl = CylStru(calculation_domain='Longitudinal Stiffened shell')
    cyl.set_material(mat_yield=355, emodule=210000, material_factor=1.15, poisson=0.3)
    cyl.set_stresses(sasd=-137.557, tQsd=78.2986, shsd=-3.8934)
    cyl.set_shell_geometry(radius=6500, thickness=19,
                           tot_length_of_shell=20000, distance_between_rings=3300)
    cyl.set_longitudinal_stiffener(hw=232, tw=12, bf=49, tf=28, spacing=680.367)
    cyl.set_panel_spacing(val=680)
    cyl.set_length_between_girder(val=3300)
    cyl.set_imperfection()
    cyl.set_fabrication_method()
    cyl.set_end_cap_pressure_included_in_stress()
    cyl.set_uls_or_als('ULS')
    cyl.set_shell_buckling_parmeters()

    for key, val in cyl.get_buckling_results().items():
        print(f'{key}: {val}')
