'''
FlatStru API tests. Reference values are taken from the original
ANYstructure implementation on identical inputs (parity verified at the
time of the code transfer).
'''
import pytest


def test_prescriptive_stiffened_panel(flat_panel_builder):
    results = flat_panel_builder('Flat plate, stiffened').get_buckling_results()

    assert results['Plate']['Plate buckling'] == pytest.approx(0.2819426684980083, rel=1e-6)
    stiffener = results['Stiffener']
    assert stiffener['Overpressure plate side'] == pytest.approx(1.4563580670007914, rel=1e-6)
    assert stiffener['Overpressure stiffener side'] == pytest.approx(1.5707510221848713, rel=1e-6)
    assert stiffener['Resistance between stiffeners'] == pytest.approx(1.3264089057658677, rel=1e-6)
    assert stiffener['Shear capacity'] == pytest.approx(0.018215841459687354, rel=1e-6)


def test_prescriptive_unstiffened_plate(flat_panel_builder):
    results = flat_panel_builder('Flat plate, unstiffened').get_buckling_results()

    assert results['Plate']['Plate buckling'] == pytest.approx(1.6482957799009763, rel=1e-6)
    assert results['Stiffener']['Overpressure plate side'] == 0


def test_prescriptive_stiffened_with_girder_runs(flat_panel_builder):
    results = flat_panel_builder('Flat plate, stiffened with girder').get_buckling_results()

    assert results['Girder']['Shear capacity'] > 0
    assert set(results) == {'Plate', 'Stiffener', 'Girder', 'Local buckling'}


def test_semi_analytical_stiffened_panel(flat_panel_builder):
    results = flat_panel_builder(
        'Flat plate, stiffened', method='SemiAnalytical S3/U3').get_buckling_results()

    assert results['method'] == 'SemiAnalytical S3/U3'
    assert results['valid prediction'] == 1
    assert results['panel family'] == 'S3'
    assert results['buckling UF raw'] == pytest.approx(1.445366740628527, rel=1e-6)
    assert results['ultimate UF raw'] == pytest.approx(0.9465671231264069, rel=1e-6)
    # UFs include the material factor of 1.15
    assert results['ultimate UF'] == pytest.approx(results['ultimate UF raw'] * 1.15, rel=1e-9)
    assert results['selected UF'] == results['ultimate UF']


def test_special_provisions(flat_panel_builder):
    results = flat_panel_builder('Flat plate, stiffened').get_special_provisions_results()

    assert set(results) == {'Plate thickness', 'Stiffener section modulus', 'Stiffener shear area'}
    for check in results.values():
        assert check['minimum'] > 0
        assert check['actual'] > 0


def test_available_methods_exclude_ml(flat_panel_builder):
    methods = flat_panel_builder().get_available_buckling_methods()

    assert 'DNV-RP-C201 - prescriptive' in methods
    assert 'SemiAnalytical S3/U3' in methods
    assert not any('ML' in method for method in methods)


def test_invalid_domain_rejected():
    from anybuckling import FlatStru

    with pytest.raises(AssertionError):
        FlatStru('Not a domain')
