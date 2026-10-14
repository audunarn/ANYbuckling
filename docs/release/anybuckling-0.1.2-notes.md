# ANYbuckling 0.1.2 — release notes

Patch release with the corrections from the calculation review
(`docs/review/CALCULATION_REVIEW.md`). The public API is backward compatible: new
parameters are optional, and no function, class or result key was removed.
Numerical results change wherever a defect was corrected. Each correction has a
regression test with an independent reference.

## DNV-RP-C201 (flat plates, stiffeners, girders)

- **Girder span:** `FlatStru.set_girder(..., span=L_G_mm)` sets it. Omitting `span` keeps the
  previous 10 m placeholder, with a warning (P1).
- **Equation fixes:**
  - (7.61)–(7.64) sniped-girder selection (P2);
  - C<sub>τe</sub> in the combined check (P5);
  - le ≤ 0.3/0.4 L<sub>G</sub> (P6);
  - p<sub>0</sub> for tension (P7);
  - τ<sub>crg</sub> for unknown L<sub>P</sub> (P8);
  - C<sub>τ</sub> (P9);
  - σ<sub>y</sub> sampling (P11);
  - (7.9) (P12);
  - c<sub>i</sub>, k<sub>p</sub> ≥ 0 and l<sub>T</sub> (P14–P16);
  - p<sub>0</sub> = 0 for tensile σ<sub>y</sub> (N2).
- **Failures reported instead of errors or wrong signs:**
  - V ≥ V<sub>Rd</sub> now fails (UF = ∞) instead of raising (P10);
  - N<sub>Sd</sub> ≥ N<sub>E</sub> now fails (UF = ∞) instead of giving a negative, passing UF (P17).
- **Plastic section modulus:** corrected for the neutral axis in the plate (P13).

## DNV-RP-C202 (cylinders, cones)

- **Hoop stress and continuity:**
  - the additional hoop stress is kept for ring-stiffened shells and cones;
  - the longitudinally stiffened shell is continuous at σ<sub>x</sub> = 0 and independent of the
    sign of σ<sub>m</sub>.
- **Column buckling:** unstiffened shells now get the column check with the 0.9/λ̄² curve.
- **Units and flags:**
  - `set_forces` uses the right units;
  - `set_end_cap_pressure_included_in_stress(is_included=False)` now does what its name says.
    The default keeps the previous results; `True` (ψ = 4) warns.
  - `set_length_between_girder` stores metres.
  - `conversions.helper_cylinder_stress_to_force_to_stress(..., length_unit=None)`.
- **Other fixes:**
  - panel σ<sub>j</sub> for curved panels;
  - L<sub>H</sub> for heavy ring frames;
  - Table 3-3 torsion 5.34;
  - fabrication flags take effect.

## Semi-analytical S3/U3

- **Ritz basis:** ⌈2a/b⌉ + 1 longitudinal half-waves, which converges both the elastic and the
  ultimate results (S1).
- **Stress gradients:** buckling drivers use RP-C201 k<sub>σ</sub>(ψ) and 6.8 (S2, U9).
- **Global stiffness:** H = D + G I<sub>t</sub>/(2s) instead of Huber √(D<sub>11</sub>D<sub>22</sub>) (S5).
- **Ultimate capacity:**
  - capped at the column-type global elastic factor (U1);
  - no lift for pressure-dominated panels, and the buckling UF is never below the ultimate UF (U2).
- **Load factors:** applied (U3).
- **Mixed-sign axial stress:** both ends are solved and the envelope is reported (U5).
- **CSR check:** uses the mm scantlings (U6).
- **Supports:** free edges are rejected; clamped and girder supports are flagged (S4).
- **Effective width:** IACS/CSR form (S7).
- **Performance:** sparse membrane contractions and field reuse. A realistic S3 panel takes about
  28 ms (19 ms in 0.1.1) on either the numba or the pure-Python path.
  `S3SolverConfig(auto_mode_count=False)` restores the old basis.

## Not included

- Open review items: P10 jump, P11 floor, N1, N3, C11, C16, C17, C19, C20, C24–C28,
  S3, S6, S8, S9, U4, U7, U8, U10, U11. See §0 and §8 of the review.
- The PULS benchmark gates have not been rerun: they need the PULS reference data, which is
  not in the repository.
