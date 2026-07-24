'''
Semi-analytical S3 (stiffened) / U3 (unstiffened) panel buckling solver.

The solver is self-contained (numpy only) and exposes dataclass inputs,
single-panel solve functions, ANYstructure-object adapters and a
benchmark / verification CLI (``python -m anybuckling.semianalytical``).
'''
from .solver import (
    S3PanelInput,
    S3Result,
    S3SolverConfig,
    U3PanelInput,
    anystructure_panel_input,
    calculate_csr_requirement,
    main,
    predict_anystructure_uf,
    predict_anystructure_uf_batch,
    predict_anystructure_uf_with_acceptance,
    row_to_s3_input,
    row_to_u3_input,
    solve_anystructure_panel,
    solve_s3_panel,
    solve_u3_panel,
    validate_s3_input,
    validate_u3_input,
)

__all__ = [
    "S3PanelInput",
    "S3Result",
    "S3SolverConfig",
    "U3PanelInput",
    "anystructure_panel_input",
    "calculate_csr_requirement",
    "main",
    "predict_anystructure_uf",
    "predict_anystructure_uf_batch",
    "predict_anystructure_uf_with_acceptance",
    "row_to_s3_input",
    "row_to_u3_input",
    "solve_anystructure_panel",
    "solve_s3_panel",
    "solve_u3_panel",
    "validate_s3_input",
    "validate_u3_input",
]
