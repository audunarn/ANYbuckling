'''
Direct tests of the semi-analytical S3/U3 solver. The pinned values are
regression values of the solver after the calculation-review fixes
(docs/review/CALCULATION_REVIEW.md); they differ from the original ANYstructure
parity values only through the aspect-ratio adaptive Ritz basis (finding S1).
With ``auto_mode_count=False`` the original parity values are reproduced
(see ``test_original_parity_with_fixed_mode_basis``).
'''
import dataclasses

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
    # Elastic buckling envelope (local plate + shear, Ritz-combined) with the
    # local plate basis m = 1..6 for L/s = 4.85 (review S1); was 1.445366740628527
    # with the fixed m <= 3 basis.
    assert result.buckling_usage_factor == pytest.approx(1.445370673941167, rel=1e-6)
    # Major-yield ultimate capacity with the m = 1..6 local plate basis (review
    # S1); the fixed m <= 3 basis gave the unconservative 0.9465671231264069.
    assert result.ultimate_usage_factor == pytest.approx(1.0415100838527238, rel=1e-6)


def test_solve_u3_panel():
    result = solve_u3_panel(row_to_u3_input(U3_ROW))

    assert result.valid is True
    assert result.invalid_reason is None
    # Elastic plate/shear buckling envelope with the plate basis m = 1..6 for
    # a/b = 4.85 (review S1); was 1.5899034146913802 with m <= 3.
    assert result.buckling_usage_factor == pytest.approx(1.5899077413352847, rel=1e-6)
    # First-yield ultimate capacity with the m = 1..6 basis (review S1); the
    # fixed m <= 3 basis gave the unconservative 0.9166210069632023.
    assert result.ultimate_usage_factor == pytest.approx(1.067753466419739, rel=1e-6)


def test_original_parity_with_fixed_mode_basis():
    # The original ANYstructure parity values are reproduced with the fixed
    # (1, 2, 3) x (1, 2) basis: the pinned changes above come from S1 alone.
    fixed = dataclasses.replace(S3SolverConfig(), auto_mode_count=False)
    s3 = solve_s3_panel(row_to_s3_input(S3_ROW), fixed)
    u3 = solve_u3_panel(row_to_u3_input(U3_ROW), fixed)

    assert s3.buckling_usage_factor == pytest.approx(1.445366740628527, rel=1e-6)
    assert s3.ultimate_usage_factor == pytest.approx(0.9465671231264069, rel=1e-6)
    assert u3.buckling_usage_factor == pytest.approx(1.5899034146913802, rel=1e-6)
    assert u3.ultimate_usage_factor == pytest.approx(0.9166210069632023, rel=1e-6)


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
