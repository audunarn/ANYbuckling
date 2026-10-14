'''
CylStru API tests. Reference values are taken from the original
ANYstructure implementation on identical inputs (parity verified at the
time of the code transfer).
'''
import pytest


def test_unstiffened_shell(cylinder_builder):
    results = cylinder_builder('Unstiffened shell').get_buckling_results()

    # The fixture uses the default set_end_cap_pressure_included_in_stress() (False: hydrostatic psi = 2,
    # DNV-RP-C202 Table 3-2 note 2). Independent hand calculation (3.4.1)-(3.4.3), (3.2.1)-(3.2.3), (3.1.3):
    # psi = 2 -> 3.3662104853009027 (psi = 4, end-cap pressure included in the axial stress -> 3.347249734791602).
    assert results['Unstiffened shell'] == pytest.approx(3.3662104853009027, rel=1e-6)
    assert results['Longitudinal stiffened shell'] is None
    assert results['Ring stiffened shell'] is None


def test_longitudinal_stiffened_shell(cylinder_builder):
    results = cylinder_builder('Longitudinal Stiffened shell').get_buckling_results()

    assert results['Unstiffened shell'] == pytest.approx(0.8557206926033175, rel=1e-6)
    assert results['Longitudinal stiffened shell'] == pytest.approx(0.9701555412318517, rel=1e-6)
    assert results['Column stability UF'] == pytest.approx(0.6967354898157304, rel=1e-6)


def test_orthogonally_stiffened_shell(cylinder_builder):
    results = cylinder_builder('Orthogonally Stiffened shell').get_buckling_results()

    assert results['Unstiffened shell'] == pytest.approx(0.8557206926033175, rel=1e-6)
    assert results['Longitudinal stiffened shell'] == pytest.approx(0.9701555412318517, rel=1e-6)
    assert results['Heavy ring frame'] == pytest.approx(0.30237783040956595, rel=1e-6)
    assert results['Column stability check'] is True
    assert results['Stiffener check']['longitudinal'] is True


def test_ring_stiffened_shell_runs(cylinder_builder):
    results = cylinder_builder('Ring Stiffened shell').get_buckling_results()

    assert results['Ring stiffened shell'] is not None
    assert results['Ring stiffened shell'] > 0
