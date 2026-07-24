'''
Direct tests of the semi-analytical S3/U3 solver. Reference values are
taken from the original ANYstructure implementation on identical inputs
(parity verified at the time of the code transfer).
'''
import pytest

from anybuckling.semianalytical import (
    S3SolverConfig,
    row_to_s3_input,
    row_to_u3_input,
    solve_s3_panel,
    solve_u3_panel,
)

S3_ROW = {
    'Length of panel': 3300.0, 'Stiffener spacing': 680.0, 'Plate thick.': 12.0,
    'Stiffener type': 'L-bulb', 'Stiffener boundary': 'Cont',
    'Stiff. Height': 260.0, 'Web thick.': 12.0, 'Flange width': 49.0,
    'Flange thick.': 28.0, 'Yield stress plate': 355.0, 'Yield stress stiffener': 355.0,
    'Axial stress': 50.0, 'Trans. stress 1': 100.0, 'Trans. stress 2': 100.0,
    'Shear stress': 5.0, 'Pressure (fixed)': 0.01, 'In-plane support': 'Integrated',
}

U3_ROW = {
    'Plate length': 3300.0, 'Plate width': 680.0, 'Plate thick.': 12.0,
    'Yield stress plate': 355.0, 'Axial stress': 50.0, 'Axial stress 2': 50.0,
    'Trans. stress': 100.0, 'Trans. Stress 2': 100.0, 'Shear stress': 5.0,
    'Pressure (fixed)': 0.01, 'In-plane support': 'Integrated',
    'Rotational support': 'SS', 'Rotational support 2': 'SS',
}


def test_solve_s3_panel():
    result = solve_s3_panel(row_to_s3_input(S3_ROW))

    assert result.valid is True
    assert result.invalid_reason is None
    assert result.buckling_usage_factor == pytest.approx(1.445366740628527, rel=1e-6)
    assert result.ultimate_usage_factor == pytest.approx(0.9465671231264069, rel=1e-6)


def test_solve_u3_panel():
    result = solve_u3_panel(row_to_u3_input(U3_ROW))

    assert result.valid is True
    assert result.invalid_reason is None
    assert result.buckling_usage_factor == pytest.approx(1.5899034146913802, rel=1e-6)
    assert result.ultimate_usage_factor == pytest.approx(0.9166210069632023, rel=1e-6)


def test_invalid_input_reported_not_raised():
    row = dict(S3_ROW, **{'Stiffener boundary': 'not-a-boundary'})
    result = solve_s3_panel(row_to_s3_input(row))

    assert result.valid is False
    assert result.invalid_reason == 'unsupported-stiffener-boundary'
    assert result.buckling_usage_factor is None


def test_result_to_dict_roundtrip():
    result = solve_s3_panel(row_to_s3_input(S3_ROW))
    as_dict = result.to_dict()

    assert as_dict['valid'] is True
    assert as_dict['buckling_usage_factor'] == result.buckling_usage_factor


def test_solver_config_is_configurable():
    config = S3SolverConfig()
    result = solve_s3_panel(row_to_s3_input(S3_ROW), config)

    assert result.valid is True
