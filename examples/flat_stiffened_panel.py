'''Flat stiffened panel checked with both calculation methods.'''
from anybuckling import FlatStru


def build_panel():
    panel = FlatStru('Flat plate, stiffened')
    panel.set_material(mat_yield=355, emodule=210000, material_factor=1.15, poisson=0.3)
    panel.set_plate_geometry(spacing=680, thickness=12, span=3300)
    panel.set_stresses(pressure=0.01, sigma_x1=50, sigma_x2=50,
                       sigma_y1=100, sigma_y2=100, tau_xy=5)
    panel.set_stiffener(hw=260, tw=12, bf=49, tf=28, stf_type='bulb', spacing=680)
    return panel


if __name__ == '__main__':
    panel = build_panel()

    panel.set_buckling_parameters(calculation_method='DNV-RP-C201 - prescriptive',
                                  buckling_acceptance='ultimate')
    print('--- DNV-RP-C201 prescriptive ---')
    for group, values in panel.get_buckling_results().items():
        print(f'{group}: {values}')

    panel.set_buckling_parameters(calculation_method='SemiAnalytical S3/U3')
    print('--- SemiAnalytical S3/U3 ---')
    results = panel.get_buckling_results()
    for key in ['buckling UF', 'ultimate UF', 'panel family', 'confidence', 'valid label']:
        print(f'{key}: {results[key]}')

    print('--- Special provisions ---')
    for check, values in panel.get_special_provisions_results().items():
        print(f'{check}: minimum {values["minimum"]:.4g}, actual {values["actual"]:.4g}')
