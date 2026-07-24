import pytest


@pytest.fixture
def flat_panel_builder():
    '''Return a builder for a reference stiffened flat panel.'''
    from anybuckling import FlatStru

    def build(domain='Flat plate, stiffened', method='DNV-RP-C201 - prescriptive'):
        panel = FlatStru(domain)
        panel.set_material(mat_yield=355, emodule=210000, material_factor=1.15, poisson=0.3)
        panel.set_plate_geometry(spacing=680, thickness=12, span=3300)
        panel.set_stresses(pressure=0.01, sigma_x1=50, sigma_x2=50,
                           sigma_y1=100, sigma_y2=100, tau_xy=5)
        if domain != 'Flat plate, unstiffened':
            panel.set_stiffener(hw=260, tw=12, bf=49, tf=28, stf_type='bulb', spacing=680)
        if domain == 'Flat plate, stiffened with girder':
            panel.set_girder(hw=500, tw=15, bf=200, tf=25, stf_type='T', spacing=700)
        panel.set_fixation_parameters()
        panel.set_buckling_parameters(calculation_method=method,
                                      buckling_acceptance='ultimate',
                                      stiffened_plate_effective_aginst_sigy=True)
        return panel

    return build


@pytest.fixture
def cylinder_builder():
    '''Return a builder for a reference cylinder.'''
    from anybuckling import CylStru

    def build(domain='Unstiffened shell'):
        cyl = CylStru(calculation_domain=domain)
        cyl.set_material(mat_yield=355, emodule=210000, material_factor=1.15, poisson=0.3)
        cyl.set_imperfection()
        cyl.set_fabrication_method()
        cyl.set_end_cap_pressure_included_in_stress()
        cyl.set_uls_or_als()
        cyl.set_length_between_girder(val=3300)
        cyl.set_panel_spacing(val=680)
        cyl.set_shell_geometry(radius=6500, thickness=19, tot_length_of_shell=20000,
                               distance_between_rings=3300)
        if 'Longitudinal' in domain or 'Orthogonally' in domain:
            cyl.set_longitudinal_stiffener(hw=232, tw=12, bf=49, tf=28, spacing=680.367)
        if 'Ring Stiffened' in domain:
            cyl.set_ring_stiffener(hw=260, tw=12, bf=49, tf=28, spacing=680)
        if 'Orthogonally' in domain:
            cyl.set_exclude_ring_stiffener()
            cyl.set_ring_girder(hw=500, tw=15, bf=200, tf=25, stf_type='T', spacing=3300)
        cyl.set_shell_buckling_parmeters()
        cyl.set_stresses(sasd=-137.557, tQsd=78.2986, shsd=-3.8934)
        return cyl

    return build
