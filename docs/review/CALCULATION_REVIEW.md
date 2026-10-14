# ANYbuckling calculation review — prescriptive and semi-analytical

| | |
|---|---|
| Reviewed revision | `a871d5a` (`main`, package version 0.1.1) |
| Date | 2026-10-14 |
| Scope | `prescriptive/plates.py` (DNV-RP-C201), `prescriptive/cylinders.py` + `conversions.py` (DNV-RP-C202), `semianalytical/solver.py` (S3/U3), and the `FlatStru` / `CylStru` API paths that feed them |
| Environment | Python 3.13, numpy 2.5.3, scipy (latest from PyPI), Linux; numba not installed (pure-Python solver paths) |
| Reproduction | `python docs/review/checks/verify_c201_plates.py`, `verify_c202_cylinders.py`, `verify_semianalytical.py` (from the repository root after `pip install -e .`) |

## 0. Status after checking the standards and fixing (update, 2026-10-14)

After the first review pass, every finding was checked against the texts of **DNV-RP-C201
(October 2010)** and **DNV-RP-C202 (January 2013)** supplied by the owner, and the confirmed
defects were fixed on this branch. Each fix has a regression test that fails without it, with
expected values from independent hand calculations, closed-form theory or properties the
standard requires (invariance, monotonicity, continuity). Each set of fixes was then reviewed
again by an independent sweep (§8). Sections 3–7 keep the original findings and the values
**before** the fixes; the table below is authoritative for their status.

| ID | Sev. | Check against the standard | Status | Commit |
|---|---|---|---|---|
| P1 | High | — (API) | Fixed: `set_girder(span=…)`; omitted → warning and the 10 m placeholder; a later call keeps an earlier span | 23068b9, 329e7dd |
| P2 | High | Confirmed, RP-C201 (7.61)–(7.64) | Fixed | 23068b9 |
| P3 | — | **Withdrawn**: 7.7.1 minimises (7.50)–(7.53) *or* (7.54)–(7.57), so a separate z\* per pressure side is allowed | No change | — |
| P4 | — | **Not a deviation**: (7.32) implies I<sub>t</sub> = ⅓(h<sub>w</sub>t<sub>w</sub>³ + b t<sub>f</sub>³) | Kept as a physics note (stocky bulb flanges) | — |
| P5 | Medium | Confirmed, (6.18)–(6.20) | Fixed | 23068b9 |
| P6 | Medium | Confirmed, 8.4.2 | Fixed | 23068b9 |
| P7 | Medium | Confirmed, (8.7) | Fixed | 23068b9 |
| P8 | Medium | (8.6) consistent; the default for unknown L<sub>P</sub> was the most favourable value | Fixed: unknown L<sub>P</sub> → τ<sub>crg</sub> = 0 (L<sub>P</sub> → ∞) | 23068b9 |
| P9 | Medium | Confirmed, (6.15) | Fixed | 23068b9 |
| P10 | Medium | 7.8 requires a reduction for V > 0.5 V<sub>Rd</sub> but gives no formula | Crash fixed (V ≥ V<sub>Rd</sub> → UF = ∞); jump at 0.5 **open** (owner decision) | 23068b9 |
| P11 | Medium | DNV-OS-C101 (not supplied) | σ<sub>y</sub> sampling fixed; 15 mm³ vs 15 cm³ floor **open** | 23068b9 |
| P12 | Low | Confirmed, (7.9) | Fixed | 23068b9 |
| P13 | Low | Theory | Fixed (and the pinned test value) | 23068b9 |
| P14 | Medium | **New**: c<sub>i</sub> = 1.0 if either stress is tensile (text to (6.18)) | Fixed | 23068b9 |
| P15 | Low | **New**: k<sub>p</sub> ≥ 0, (6.10) | Fixed; zero resistance fails (UF = ∞) | 23068b9, 329e7dd |
| P16 | Low | **New**: l<sub>T</sub> = min(0.4 l, support distance) etc., (7.68)/(7.69), (8.12)/(8.13) | Fixed | 23068b9 |
| N1 | Low | **New (sweep)**: C<sub>yG</sub> for tensile σ<sub>y</sub> "may be taken as 1.0" (8.21); code uses a smaller value | **Open** (conservative) | — |
| N2 | Low | **New (sweep)**: p<sub>0</sub> = 0 when σ<sub>y</sub> is nowhere compressive, (7.9)/(7.10) | Fixed | 329e7dd |
| N3 | Low | **New (sweep)**: 'Local buckling' returns allowable dimensions (9.1), not UFs; flat-bar web not checked against (9.2) | Docstring fixed; flat-bar web **open** | 329e7dd |
| C1 | High | Internal | Fixed | 46cd5a6 |
| C2 | High | Internal (Sec 4.2.3) | Fixed; stress/force input recorded explicitly by the API | 46cd5a6, d301655 |
| C3 | High | (3.2.2), (3.6.4), (3.6.7) | Fixed: term-by-term λ<sub>s</sub>², C → ρξ as s<sub>e</sub> → 0 | 46cd5a6 |
| C4 | High | Confirmed, (2.2.3), (3.2.5) | Fixed (geometries 3 and 7) | 46cd5a6 |
| C5 | High | Confirmed, 3.4.1 b), 3.8 | Fixed: column check (3.8.1)–(3.8.12) for geometries 1 and 2 | 46cd5a6 |
| C6 | High | Confirmed, (3.8.6) | Fixed | 46cd5a6 |
| C7 | High | Internal (2.2.2)–(2.2.7) | Fixed | 46cd5a6 |
| C8 | Medium (was High) | 3.5.2.7 is an allowed alternative to I<sub>h</sub> | Fixed (index) | 46cd5a6 |
| C9 | Medium | (3.1.3) is defined for γ<sub>M0</sub> = 1.15 only | Fixed (whole-curve scaling, continuous) | 46cd5a6 |
| C10 | Medium | Confirmed, (3.10.5) | Fixed | 46cd5a6 |
| C11 | Medium | (3.10.1)/(3.10.2) are not pass/fail: failing them means f<sub>r</sub> = f<sub>T</sub> | **Open**: needs redesign, not re-enabling | — |
| C12 | Medium | Confirmed, 3.3.2 | Fixed | 46cd5a6 |
| C13 | Medium | Confirmed, Table 3-2 notes 1/2 | Fixed: `is_included=True` → ψ = 4, default `False` → ψ = 2 (previous default results kept) | 46cd5a6, d301655 |
| C14 | Medium | Confirmed, 3.5.2.6 | Fixed | 46cd5a6 |
| C15 | — | **Withdrawn**: 3.8.2 defines f<sub>ak</sub> with the actual σ<sub>m</sub>, σ<sub>h</sub>, τ | No change | — |
| C16, C17, C19, C20 | Low | — | **Open** (C20's σ<sub>h</sub> double count at line 692 fixed) | 46cd5a6 |
| C18 | Low | Confirmed, (3.5.4) | Fixed | 46cd5a6 |
| S1 | High | Theory | Fixed: `auto_mode_count` (≈ a/b half-waves) | c1f42ef |
| S2 | High | RP-C201 6.6, 6.8 | Fixed (interim): equivalent drivers from k<sub>σ</sub>(ψ) and 6.8, ψ < −2 bounded by k<sub>σ</sub>(−2) | c1f42ef, d029b55 |
| S3, S5, S6, S8, S9 | Med./Low | Theory / ECOSYSTEM_GUIDE | **Open** (empirical constants and stiffness choices) | — |
| S4 | Medium | Internal | Fixed: free edge rejected; CL and GL/GT flagged in diagnostics | c1f42ef |
| S7 | Medium | IACS CSR | Fixed | c1f42ef |
| U1 | High | Theory | Fixed: capped at the column-type global elastic factor + stability-loss check; confidence "low" | c1f42ef |
| U2 | High | ECOSYSTEM_GUIDE | Fixed: lift removed | c1f42ef |
| U3 | High | Internal | Fixed: load factors applied; pressure side recorded only (plate side modelled) | c1f42ef |
| U4 | High | Theory | **Open** (sniped stiffeners under lateral pressure) | — |
| U5 | High | Internal | Fixed: both ends solved, envelope reported | c1f42ef |
| U6 | High | Internal | Fixed | c1f42ef |
| U7, U8, U10, U11 | Med./Low | Theory / ECOSYSTEM_GUIDE | **Open** | — |
| U9 | Medium | Internal (equilibrium) | Fixed | c1f42ef |

**Behaviour changes users will notice**

- `FlatStru.set_girder(..., span=L_G_mm)`: girder checks now use the real girder span; omitting
  it warns.
- `CylStru.set_end_cap_pressure_included_in_stress(is_included=False)` now means what its name
  says, and its default is `False`, which keeps the previous default results.
- Unstiffened cylinders and curved panels now report `Column stability check` /
  `Column stability UF`.
- SemiAnalytical results change for panels with L/s > 3 (S1), stress gradients (S2, U9), load
  factors ≠ 1 (U3), mixed-sign axial stress (U5), pressure-dominated panels (U2) and slender
  columns (U1). Example: README S3 panel, ultimate UF raw 0.947 → 1.042; U3 test row
  0.917 → 1.068. Solves are slower for long panels (S3 up to ≈ 0.1 s at L/s = 12, previously
  ≈ 0.01 s).

## 1. How the review was done, and its limits

* Every formula path listed in the scope was read line by line and compared with the
  standard. Key intermediate values were recomputed by an independent hand
  implementation, and each suspected defect was reproduced numerically through the
  public API before being listed here. Suspected issues that did not reproduce were
  dropped (for example, a suspected sign flip of the C201 interaction equations when
  N<sub>Sd</sub> > N<sub>E</sub> turned out to be harmless because N<sub>Sd</sub>/N<sub>ks,Rd</sub> already exceeds 1).
* The first pass was made without the texts of the standards (network access to DNV was
  blocked); the owner then supplied DNV-RP-C201 (Oct 2010) and DNV-RP-C202 (Jan 2013) and every
  finding was re-checked against them (§0). Each finding states its basis:
  * **Internal** – the code contradicts itself (wrong variable, dead variable, inverted
    branch, unit slip, a setting that is silently ignored). These do not depend on how
    the standard is worded.
  * **Standard** – the code deviates from the standard (now checked against the text; see §0).
  * **Theory** – the semi-analytical solver deviates from classical plate/beam theory
    (closed-form results that are not in doubt).
* Severity: **High** = can report a passing UF where the correct UF fails, or the number
  is plainly wrong in normal use. **Medium** = error confined to a sub-domain, or
  conservative but material, or an input that is silently ignored. **Low** = robustness,
  dead code, cosmetic.
* Evidence: every number in this report is reproduced by the three committed scripts,
  except numbers marked **†**, which were measured with the reviewers' larger exploratory
  scripts (not committed) and checked only at source level. Findings marked "source" were
  established by reading the code.
* The existing test suite (104 tests, all passing) checks **parity with the original
  ANYstructure implementation**, not correctness. No test compares with an independent
  reference, so none of the defects below is caught; one test (`test_plastic_sec_mod`)
  pins a value that is wrong by a factor of 19 (P13).

## 2. Summary

| Module | Verified correct | High | Medium | Low |
|---|---|---|---|---|
| DNV-RP-C201 plates / stiffeners / girders | Most formula chains (§3.1) | 4 | 7 | 2 |
| DNV-RP-C202 cylinders / cones | Tables 3-1…3-3, ring stresses, 3.9, cones (§4.1) | 8 | 6 | 6 |
| Semi-analytical S3/U3 – elastic core | Integrals, membrane field, sections, Euler column (§5.1) | 2 | 6 | 6 |
| Semi-analytical S3/U3 – ultimate path, adapters, API | Newton solve, stress recovery, pressure limits, tension/yield limits (§6.1) | 6 | 4 | 6 |

**Overall verdict.** The prescriptive DNV-RP-C201 and DNV-RP-C202 implementations are largely
faithful to the standards along their main paths — the hand calculations reproduce the
README results exactly — but each has a handful of defects that change pass/fail in ordinary
use. The semi-analytical S3/U3 solver has correct numerical building blocks (integrals,
membrane field, Newton solve), but its ultimate-capacity logic is **not fit for acceptance
decisions in its current form**: it can report UF 0.33 for a stiffener that fails Euler
buckling by a factor of five (U1).

Most urgent, because they change pass/fail in ordinary use:

1. **U1** – SemiAnalytical ultimate capacity can exceed the global (Euler) buckling load: a 7 m flat-bar stiffener with Euler stress 20.5 MPa under 100 MPa reports selected UF 0.33 (RP-C201: 6.52).
2. **U2 / U3 / U5** – SemiAnalytical: more lateral pressure can *lower* the ultimate UF to below σ/f<sub>y</sub>; load factors are ignored; a 300 MPa tension end is ignored (UF 0.12).
3. **P1** – `FlatStru` fixes the girder span at 10 m; girder UFs do not respond to the real girder.
4. **P2** – sniped girders under transverse compression: wrong interaction equations (UF 0.17 instead of 1.14).
5. **C1 / C2** – the user's additional hoop stress is discarded for ring-stiffened shells and for cones.
6. **C3 / C4** – longitudinally stiffened shells collapse to a near-zero slenderness when σ<sub>x</sub> = 0, and depend on the sign of the bending stress.
7. **C5 / C6** – no column-buckling check for unstiffened tubes; column curve 11 % unconservative beyond λ̄ = 1.34.
8. **C7** – `CylStru.set_forces` produces bending and torsion stresses 1000× too large.
9. **S1 / S2** – the semi-analytical solver's fixed 3×2 Ritz basis and averaged stress gradients make the elastic buckling factor unconservative for L/s > 3 (+28 % at L/s = 5, +352 % at 12) and blind to in-plane bending.

## 3. DNV-RP-C201 — `prescriptive/plates.py` and `FlatStru`

### 3.1 Verified

An independent hand calculation of the README stiffened panel (bulb 260×12 / 49×28,
s = 680, t = 12, l = 3300, σ<sub>x</sub> = 50, σ<sub>y</sub> = 100, τ = 5, p = 0.01 MPa) reproduces the code
exactly for λ<sub>p</sub>, C<sub>xs</sub>, σ<sub>y,R</sub>, C<sub>ys</sub>, I<sub>s</sub>, I<sub>po</sub> and I<sub>z</sub>
(`verify_c201_plates.py`, block V0). Note that in that example σ<sub>y</sub> = 100 MPa exceeds
σ<sub>y,R</sub> = 86.7 MPa, so the effective plate width is zero; the README example is an
overloaded panel.

Checked and consistent with RP-C201 (2010): Sec 5 lateral-pressure plate check (ψ<sub>x</sub>, ψ<sub>y</sub>);
6.2/6.6 longitudinal compression incl. k<sub>σ</sub>(ψ) and C<sub>x</sub>; 6.3/6.7 transverse compression
(λ<sub>c</sub>, μ, κ, k<sub>p</sub>, h<sub>α</sub>, the l<sub>1</sub> = min(l/4, s/2) rule and the 0.75σ<sub>max</sub> floor);
6.4 k<sub>l</sub>, λ<sub>w</sub>; 6.5 interaction form and c<sub>i</sub>; 7.2 N<sub>Sd</sub>, k<sub>c</sub>, m<sub>c</sub>, C<sub>0</sub>, p<sub>0</sub>,
τ<sub>crl</sub>, τ<sub>crg</sub>, τ<sub>tf</sub>; 7.3 effective width (compression and tension forms); 7.4 k<sub>sp</sub>; 7.5.1
column curve, f<sub>E</sub>, i<sub>e</sub>, μ for plate and stiffener side; 7.5.2 f<sub>ET</sub> structure, β, C, η,
f<sub>Ep</sub>, λ<sub>e</sub>, f<sub>Epx/y/τ</sub>, c, I<sub>po</sub>, I<sub>z</sub>, h<sub>s</sub>; 7.6 I<sub>p</sub>, τ<sub>crs</sub>, τ<sub>Rd</sub>, u; 7.7.1/7.7.2
equation forms (7.50–7.63), M<sub>1</sub> = ql²/12, M<sub>2</sub> = ql²/24, l<sub>k</sub> and p<sub>f</sub>; 7.7.3 resistance
parameters; 8.5 L<sub>GT0</sub> and f<sub>ET,G</sub>; elastic section properties (I<sub>y</sub>, W, centroid).

### 3.2 Findings

**P1 — Girder span L<sub>G</sub> is fixed at 10 m in `FlatStru`; `set_girder(spacing=…)` does nothing.** High · Internal
`api.py:264, 336, 376` set `girder_lg = 10` ("placeholder value") and no public setter exists;
`plates.py:1901, 2268` use it as L<sub>G</sub> (also in k<sub>g</sub> / τ<sub>crg</sub> for the stiffener check).
`Girder.spacing` is never read by `girder()`. Girder UFs for the README panel with a girder:
0.18 (L<sub>G</sub> = 4 m), 1.90 (10 m), 13.7 (16 m) — but through the API always 1.90, whatever is passed.

**P2 — Sniped girders: equation-set selection for pressure on the girder side is inverted.** High · Internal
`plates.py:2481–2486`. When q·L<sub>G</sub>²/8 < N<sub>y,Sd</sub>·z\*, RP-C201 uses (7.62)/(7.63); the code uses
(7.60)/(7.61), whose bending term is then negative. (The stiffener version at
`plates.py:2174–2178` is correct.) σ<sub>y</sub> = 40, no pressure: code UF 0.175, correct 1.144.

**P3 — z\* is optimised separately for the plate-side and stiffener-side equation sets.** ~~High · Standard~~ **Withdrawn** after checking the text: RP-C201 7.7.1 selects z\* to minimise (7.50)–(7.53) *or* (7.54)–(7.57).
`plates.py:2201, 2220` (stiffeners) and `2506, 2520` (girders). RP-C201 7.7.1 selects one z\* so that
the maximum of (7.50)–(7.57) is minimised; z\* is a property of the end connection and cannot
differ between the two pressure directions. With the default `pressure_side='both sides'`,
reporting min<sub>z</sub> A(z) and min<sub>z</sub> B(z) instead of min<sub>z</sub> max(A, B) is unconservative:
README panel 1.571 → 1.710; flat bar 200×14, σ<sub>x</sub> = 120, p = 0.1: 0.878 → 1.148.

**P4 — St. Venant torsion constant omits the (1 − 0.63 t<sub>f</sub>/b<sub>f</sub>) flange factor.** **Not an RP deviation** after checking the text: (7.32) implies the thin-wall I<sub>t</sub> used by the code. Kept as a physics note for stocky bulb idealisations.
`plates.py:651`: I<sub>t</sub> = (h<sub>w</sub>t<sub>w</sub>³ + b<sub>f</sub>t<sub>f</sub>³)/3. RP-C201 7.5.2:
I<sub>t</sub> = ⅓(h<sub>w</sub>t<sub>w</sub>³ + b<sub>f</sub>t<sub>f</sub>³(1 − 0.63 t<sub>f</sub>/b<sub>f</sub>)). For the bulb-as-angle idealisation
(49×28, t<sub>f</sub>/b<sub>f</sub> = 0.57) I<sub>t</sub> is 34 % too high and f<sub>T</sub> up to 20 % too high. README bulb,
σ<sub>y</sub> = 60: stiffener-side UF 0.524 → 0.606 (+16 %). T 400×12/150×20: +0.4 %. Also used by the
cylinder torsional-buckling check (`cylinders.py:633`).

**P5 — Sec 6.5 combined check never uses C<sub>τe</sub>.** Medium (unconservative) · Internal
`plates.py:1861–1862`: `tauRd_comb = Ctaue*…` is immediately overwritten by
`tauRd_comb = tauRd if sysd>0 else tauRd` (both branches identical). For compressive σ<sub>y</sub> the
elastic C<sub>τe</sub> (no post-buckling shear) must be used. Unstiffened plate s/t = 100, λ<sub>w</sub> = 1.36:
UF 0.273 → 0.330 (+21 %). Affects the unstiffened-plate domain only.

**P6 — Girder effective width: the shear-lag cap l<sub>e</sub> ≤ 0.3 L<sub>G</sub> (0.4 L<sub>G</sub> sniped) is computed but never used.** Medium · Internal
`plates.py:2327–2339`: `eff_width_sec_mod` is assigned and never read; W<sub>ep</sub>, W<sub>eG</sub> use the
uncapped method-1 width. L<sub>G</sub> = 6 m: W<sub>ep</sub> overestimated by 74 %.

**P7 — Girder p<sub>0</sub> for tensile σ<sub>x</sub> has two typos.** Medium · Internal
`plates.py:2397`: `(l - s/Lg)` instead of `(1 - s/Lg)` and `gammaM` instead of `fy`, versus the
compression form two lines above. The shear-induced lateral load on the girder is ~10⁶× too small,
i.e. ignored, whenever σ<sub>x</sub> is tensile.

**P8 — τ<sub>crg</sub> defaults to 0.6 f<sub>y</sub> when the panel length L<sub>p</sub> is not given.** Medium · Standard
`plates.py:2299–2304`: L<sub>p</sub> = None → α<sub>τ1</sub> = 0 → τ<sub>crg</sub> = 0.6 f<sub>y</sub>, the most favourable value
(L<sub>p</sub> → ∞ should give τ<sub>crg</sub> → 0). This removes the shear contribution C·τ to the girder load p<sub>0</sub>.
Girder UF 2.07 with L<sub>p</sub> unset vs 2.38 with L<sub>p</sub> = 12 spans.

**P9 — C<sub>τ</sub> middle branch.** Medium (conservative) · Standard
`plates.py:1842–1843`: 1 − 0.675(λ<sub>w</sub> − 0.8) for λ<sub>w</sub> < 1.25. RP-C201 6.4 (ENV 1993 post-critical
form) is 1 − 0.625(λ<sub>w</sub> − 0.8) for λ<sub>w</sub> ≤ 1.2, which is continuous with 0.9/λ<sub>w</sub>; the coded form
jumps at 1.25 (0.696 → 0.720). Up to 2.9 % conservative.

**P10 — Reduced web thickness for V<sub>Sd</sub>/V<sub>Rd</sub> ≥ 0.5.** Medium (robustness) · Standard (formula) / Internal (crash)
`plates.py:2042, 2351`: t<sub>w,red</sub> = t<sub>w</sub>(1 − V/V<sub>Rd</sub>) halves the web the moment V/V<sub>Rd</sub> reaches 0.5
(UF 5.58 → 15.2 across V/V<sub>Rd</sub> = 0.493 → 0.501) and reaches zero at V/V<sub>Rd</sub> = 1, where
`get_buckling_results()` raises `ZeroDivisionError`. Please confirm the RP formula; any
continuous form (e.g. one that equals t<sub>w</sub> at 0.5) removes the jump, and V/V<sub>Rd</sub> ≥ 1 should
return a failing UF rather than raise.

**P11 — Special provisions (DNV-OS-C101).** Medium · Internal
`plates.py:1101–1102, 1139–1140`: σ<sub>y</sub> is interpolated from σ<sub>y2</sub> towards σ<sub>y1</sub> by l<sub>1</sub>/l, i.e. sampled near
the σ<sub>y2</sub> end. The same stress field gives Z<sub>min</sub> = 442 043 or 760 087 mm³ depending on which end
is labelled 1. (`get_plate_efficent_b`, line 1271, does this correctly.) The minimum-section
floor `15**3/1000**3` m³ is 3 375 mm³; the OS-C101 floor is 15 × 10³ mm³.

**P12 — p<sub>0</sub> uses ψ = 1 regardless of the σ<sub>y</sub> gradient.** Low (conservative) · Internal
`plates.py:1915` hard-codes `shear_ratio_trans = 1`, so (0.6 + 0.4ψ) is always 1.0.

**P13 — `Structure.get_plasic_section_modulus` is wrong when the plastic neutral axis lies in the plate.** Low (not used by the buckling path) · Theory
`plates.py:501–511` assumes the neutral axis is in the web and mixes reference faces. Test panel
`obj_dict`: code 1.78 × 10⁻² m³, hand 9.41 × 10⁻⁴ m³ (19×). `tests/test_structure.py::test_plastic_sec_mod`
pins the wrong value.

## 4. DNV-RP-C202 — `prescriptive/cylinders.py`, `conversions.py` and `CylStru`

### 4.1 Verified

The cylinder reviewer's independent re-implementation† (shell and panel buckling, ring and stiffener checks) reproduces the
README/test results exactly: curved-panel UF 0.85550, longitudinally stiffened UF 0.96997, heavy
ring frame UF 0.30238. Checked and consistent: 3.1 γ<sub>M</sub> for γ<sub>M0</sub> = 1.15 and 1.10; 3.2.1–3.2.6
(σ<sub>j</sub>, compressive parts, λ<sub>s</sub>, f<sub>ks</sub>); Tables 3-1, 3-2, 3-3 (ψ, ξ, ρ, C, both ρ variants);
Z<sub>s</sub>, Z<sub>l</sub>; (3.4.4), (3.4.5) and their trigger limits; ring stress distribution (α, β, l<sub>eo</sub>, ζ,
σ<sub>h</sub>, σ<sub>h,R</sub>); A<sub>req</sub>, I<sub>x</sub>, I<sub>xh</sub>, I<sub>h</sub>, r<sub>0</sub>, I<sub>y</sub>; 3.9 f<sub>ET</sub>, β, C, μ, f<sub>T</sub>, l<sub>T</sub> = π√(r h<sub>w</sub>);
stocky-profile criteria; 3.8 Euler stress, i<sub>c</sub>, need-to-check test; cone r<sub>e</sub>, l<sub>e</sub>, t<sub>e</sub> and the
stress transformations; MPa/Pa and mm/m chains for stress input; γ<sub>M</sub> applied exactly once per
UF. Limiting cases: r → ∞ recovers the flat-plate f<sub>Ea</sub>; long-cylinder axial
f<sub>Ea</sub> = ρ·0.605·E t/r; pure hoop pressure (short and long, 0.25E(t/r)²); long-cylinder torsion.

### 4.2 Findings

**C1 — Ring-stiffened shells (geometries 5/6) drop the additional hoop stress σ<sub>h,Sd</sub>.** High · Internal
`cylinders.py:707–709` builds the ring-modified hoop stress without `+ self._shsd`, and
`unstiffened_shell` reads that entry for ring-stiffened geometries (`869–872`). σ<sub>a</sub> = −40,
τ = 10, σ<sub>h,add</sub> = −50: ring-stiffened shell UF 0.738 (same as with no hoop stress); the same
shell as geometry 1 gives 2.624.

**C2 — Cones with stress input discard the user's σ<sub>h,Sd</sub>.** High · Internal
`cylinders.py:809–810` sets `_psd = 0` and replaces `_shsd` by the pressure-only hoop stress
from `conical_stress_state`. UF 0.7329 for σ<sub>h</sub> = 0, −50 and −100 MPa alike.

**C3 — Longitudinally stiffened shell (3.6) collapses when σ<sub>x</sub> = 0.** High · Internal
`cylinders.py:1371`: s<sub>e</sub> ∝ |σ<sub>x</sub>| = 0 → C = 0, f<sub>Ea</sub> = 0 → the guard at `1445–1446` sets
λ<sub>s</sub>² = 0 and the hoop and shear terms vanish. τ = 50, σ<sub>h</sub> = −100: UF 0.445 at σ<sub>x</sub> = 0,
1.483 at σ<sub>x</sub> = −0.5 MPa.

**C4 — Longitudinally stiffened full shell depends on the sign of σ<sub>m</sub>.** High · Theory
`cylinders.py:1433–1452` treats only negative σ<sub>m</sub> as compressive, but a full cylinder in bending
always has a compressed side. σ<sub>a</sub> = −20, σ<sub>h</sub> = −100: UF 1.029 with σ<sub>m</sub> = +50, 1.318 with −50.

**C5 — No column-buckling check (3.8) for unstiffened shells.** High · Standard
`cylinders.py:478–579` calls `column_buckling` only for geometries 3–8; called directly it
raises `ZeroDivisionError` (`1632`). Tube r = 500, t = 12, L = 40 m, σ<sub>a</sub> = 150: (L/i<sub>c</sub>)² = 12 800 ≫
2.5E/f<sub>y</sub> = 1 479, Euler stress 162 MPa; the code reports shell UF 0.506 and no column result
(reviewer's estimate of the column UF: 1.19†).

**C6 — Column curve for λ̄ > 1.34 lacks the 0.9 factor.** High · Standard
`cylinders.py:1512–1514`: f<sub>kc</sub> = f<sub>ak</sub>/λ̄²; RP-C202 3.8.2 gives 0.9 f<sub>ak</sub>/λ̄², continuous with
(1 − 0.28λ̄²) at 1.34. The coded curve jumps from 0.497 to 0.557 at λ̄ = 1.34 and is 11 %
unconservative beyond it.

**C7 — `CylStru.set_forces`: σ<sub>m</sub> and τ<sub>T</sub> are 1000× too large.** High (conservative but unusable) · Internal
`conversions.py:36–38` (and the inverse at 60–61) use factors for mm geometry with r, t in metres.
M = T = 10 000 kNm, r = 6.5 m, t = 19 mm: σ<sub>m</sub> = 3 971 MPa (hand 3.965), τ<sub>T</sub> = 1 983 MPa
(hand 1.983). Axial force and shear force are converted correctly.

**C8 — Ring stiffener I<sub>h</sub> requirement silently excluded for geometries 5/6.** Medium · Internal (3.5.2.7 is an allowed alternative to I<sub>h</sub>)
`cylinders.py:676–678` appends `nan` for the absent longitudinal stiffener, shifting `shRsd` so
`ring_stiffened_shell` (`1257–1266`) reads `nan` for the ring; `Ih > 0` is then False and I<sub>h</sub> is
dropped from I<sub>R</sub>. Reviewer's example† (p = −0.3, σ<sub>a</sub> = −20, τ = 5): ring UF 3.60 vs 22.1 with the
correct index; the index shift itself is confirmed at source level. (If the refined 3.5.2.7 check is accepted as an alternative to I<sub>h</sub>, 3.60 may be
admissible, but it is inconsistent with geometry 7, which includes I<sub>h</sub>.)

**C9 — γ<sub>M</sub> for material factors other than 1.15/1.10.** Medium · Internal (source)
`cylinders.py:961, 1096, 1474`: `0.85+0.6*lambda_s*(mf/1.15)` scales only the 0.6λ<sub>s</sub> term, while
the λ<sub>s</sub> > 1 branch scales the whole 1.45. γ<sub>M</sub> is discontinuous at λ<sub>s</sub> = 1 (and non-monotonic
for mf < 1.15).

**C10 — Ring flange-eccentricity check compares e<sub>f</sub> in metres with a dimensionless limit.** Medium · Internal (source)
`cylinders.py:1765` uses `get_flange_eccentricity()` (metres, not divided by t<sub>w</sub>), so the check always
passes.

**C11 — Ring-frame local-buckling checks are disabled.** Medium · Internal (source)
`cylinders.py:1802–1804` returns `True` unconditionally (original checks commented out).

**C12 — Unstiffened curved panels ignore the panel width.** Medium (conservative) · Internal
`cylinders.py:856–857` uses s = min(l, 2πr) for geometries 2/6; `set_panel_spacing` only selects the
reported UF. UF 3.1737 for s = 680, 1000, 2000 mm (reference 0.855, 1.292, 2.625†).

**C13 — End-cap flag has no effect.** Medium (conservative) · Internal
`cylinders.py:1061` compares with the string `'not included in axial stresses'`; the API stores a
bool, so ψ for hydrostatic pressure is always used.

**C14 — Fabrication method has no effect.** Medium (conservative) · Internal
`cylinders.py:1255–1256` compares with `1`; the API stores `'Fabricated'` / `'Cold formed'`, so
f<sub>r</sub> = 0.9 f<sub>T</sub> is always used.

**C15–C20 — Low (†/source).** ~~f<sub>ak</sub> for the column check includes σ<sub>m</sub>, which the interaction formula counts
again~~ (C15 **withdrawn**: RP-C202 3.8.2 defines f<sub>ak</sub> with the actual σ<sub>m</sub>, σ<sub>h</sub>, τ applied). Column section properties use the
dimensions of the last stiffener object looped over (`1575–1578`). η in 3.9 mixes panel f<sub>ks</sub> with
cylinder γ<sub>M</sub> (`1563–1573`, negligible). Shell effective width l<sub>e</sub> not capped at l (`634/639`).
A zero f<sub>ak</sub> is reported as "Column stability UF = 0". G = 80 769.2 and ν = 0.3 hard-coded.

## 5. Semi-analytical S3/U3 — elastic / Ritz core (`solver.py` ≈ 1371–4160)

### 5.1 Verified

(All by the elastic-core reviewer's scripts†, plus S1/S2 here.) Section area, neutral axis and I<sub>x</sub> for flat bar, T, angle and bulb (exact); I<sub>t</sub>; D = Et³/12(1 − ν²);
D<sub>x</sub> = EI/s; every closed-form Ritz integral (sin·cos, half-range, shear coupling incl. the
m + m′ odd / n + n′ odd selection rule, pressure work 4abp/(mnπ²)) against quadrature to ~10⁻¹⁵;
membrane field (compatibility, Airy stresses, energy Gram matrix, force and tangent kernels) to
10⁻¹¹ by finite differences; plate eigenvalues within the basis for uniaxial, transverse and
biaxial loading (exact); Euler column π²EI/(AL²) for all profile types (exact). Units are
N, mm, MPa; compression positive; E = 210 000, ν = 0.3.

### 5.2 Findings

**S1 — The Ritz basis is fixed at m ≤ 3, n ≤ 2.** High · Theory
`solver.py:200–201`, used at `2762–2763` (S3) and `2896–2897` (U3). Axial plate buckling needs
m ≈ a/b half-waves. U3 elastic factor vs the classical SS plate: exact to a/b = 3, then +8.5 %
(4), +28 % (5), +131 % (8), +352 % (12; inside the accepted domain limit `max_aspect_ratio = 12`).
The ANYstructure adapters pass span × spacing, so a/b ≈ 3–8 is typical. On the README S3 panel
the raw ultimate UF rises from 0.947 to 1.078 with m ≤ 7, n ≤ 3. The same truncation affects
wide U3 plates under σ<sub>y</sub> and σ<sub>x</sub> + τ.

**S2 — Stress gradients are replaced by their mean in the elastic analysis.** High · Theory
`normalized_load_components` (`1384–1385`) uses `axial_stress` / `mean_transverse_stress`
= ½(σ<sub>1</sub> + σ<sub>2</sub>). Pure in-plane bending (ψ = −1, classical k ≈ 23.9) gives no buckling driver at all;
ψ = −0.5 is 18 %† unconservative; σ<sub>y</sub> varying along an S3 stiffener (ψ = 0) 26–41 %†.

**S3 — `flanged_local_plate_restraint_factor = 1.10`.** Medium-High · Theory / ECOSYSTEM_GUIDE
`solver.py:214`, applied at `1549–1556` to D<sub>11</sub>, D<sub>22</sub> of every local plate mode for flanged
stiffeners. An uncalibrated +10 % stiffness that does not vanish as the stiffener approaches
tripping. README S3 panel (with a converged basis): ultimate UF 1.078 → 1.149 when set to 1.0.

**S4 — Support conditions accepted but ignored.** Medium · Internal
U3 rotational supports `CL`/`FS` pass validation (`2100–2103`) but never enter the mechanics: SS and
CL give identical results. CL treated as SS is conservative; a free edge (FS) treated as SS is
unconservative (reviewer: elastic UF 0.571 instead of ≈ 5.4†). Likewise the in-plane support type
(Integrated / GL / GT) only feeds `_support_membrane_factor`, which is unused when the membrane
field exists, so all three give identical UFs†.

**S5 — Huber cross-rigidity √(D<sub>11</sub>D<sub>22</sub>) for the stiffened (global) family.** Medium · Theory
`solver.py:1565–1566`. For an open-stiffener plate H ≈ D + GI<sub>t</sub>/(2s); the Huber form overstates the
transverse global candidate by 76–94 %† (did not govern in tested cases). "Invented stiffness"
in ECOSYSTEM_GUIDE terms.

**S6 — Tripping formula is empirical.** Medium · Theory
`_stiffener_torsional_buckling` (`3690–3745`): 0.65·(G I<sub>t</sub> + E I<sub>w</sub>(π/0.7L)²)/I<sub>p</sub>. No plate
rotational restraint, no minimisation over half-waves, angle/bulb flange eccentricity ignored in
I<sub>w</sub>, and the 0.7L effective length for continuous stiffeners has no basis (continuous stiffeners in
uniform compression behave as simply supported with alternating bays). Against RP-C201 7.5.2 with
β = 1 the result is 12–25 %† unconservative for continuous T-bars and conservative elsewhere.

**S7 — `effective_stiffener_plate_width` exponent has the wrong sign.** Medium (latent) · Standard
`solver.py:1466`: 1.12/(1 + 1.75 r<sup>1.6</sup>); IACS/CSR: 1.12/(1 + 1.75 r<sup>−1.6</sup>). At l<sub>eff</sub>/s = 2.89 the
coefficient is 0.106 instead of 0.848. Off by default (`use_effective_stiffener_width=False`).

**S8 — S3 shear knock-downs.** Medium (conservative) · ECOSYSTEM_GUIDE
`s3_shear_buckling_capacity_factor = 0.75` (`3360`, `3601`), plate τ applied as web in-plane shear
(`3595–3602`) and a `(u_p^1.2 + u_w^1.2)^(1/1.2)` interaction stack up to 36–40 %† below the SS plate
in pure shear.

**S9 — Low (†/source).** Shear-only imperfection falls back to mode (1,1), which can be orthogonal to the critical
shear mode (`2708–2721`). Absolute `EPS = 1e-12` floors on dimensional quantities drop membrane
harmonics for L ≳ 3.5 m (`2274`, `2375`, `2411`; < 10⁻⁵ effect on UF in tests). Flat-bar web check
uses k = 4 instead of ≈ 0.425 (`3584`; masked by tripping). A pressure mode with ×5 stiffness
(`224`, `2788–2793`) and a global n ≥ 1 family with EI/s on shapes that vanish at the stiffeners
(`1558–1561`) — no measurable effect, but invented stiffness. Dead fallbacks
`nonlinear_membrane_factor = 0.75` and `_support_membrane_factor = 0.9`.

### 5.3 Empirical constants (ECOSYSTEM_GUIDE: "hidden corrections, empirical gate-tuning, invented stiffness … are prohibited")

| Constant | Where | Effect | Basis |
|---|---|---|---|
| 1.10 local plate restraint | 214, 1549–1556 | +10 %, unconservative | none |
| 0.75 S3 shear capacity | 215, 3360, 3601 | −25 %, conservative | none |
| 1.20 plate/web interaction exponent | 213, 3667–3672 | conservative | none |
| 0.65 torsional restraint, 0.70 L tripping length | 220, 3730–3732 | mixed | none |
| 1.2 sniped eccentricity | 217, 4590 | web demand 2.5× nominal | partly |
| 5.0 pressure-mode stiffness | 224, 2788–2793 | no observed effect | none |
| b/200, L/1000 imperfections; p s⁴/384D | 208–209, 2738 | — | DNV-CG-0128 Sec 6; clamped strip — justified |

## 6. Semi-analytical S3/U3 — ultimate path, adapters and API

### 6.1 Verified

* README trace `FlatStru` → `S3PanelInput`: geometry, E, f<sub>y</sub>, stresses and pressure arrive with the
  right units (MPa, mm; pressure MPa → kPa → MPa) and sign (compression positive on both sides);
  `api.py:195–196` applies γ<sub>M</sub> exactly once.
* Newton solve†: tangent equals the central-difference Jacobian of the residual to 2 × 10⁻¹¹ and is
  symmetric; the residual is the exact gradient of the membrane energy; per-component step limiting;
  non-convergence after cut-backs is returned as an explicit `non-convergence` result.
* Stress recovery†: membrane stresses match the single-mode closed form (±56.03 MPa); von Mises
  formula correct; the 9 × 9 critical-point grid misses the 161 × 161 maximum by ≤ 0.12 %.
* Pressure limits†: stiffener preload p s L²/12/(Z<sub>top</sub> f<sub>y</sub>), p<sub>Fs</sub>, web shear and plate strip all
  match beam theory.
* Limiting cases†: tension (U3/S3, uniaxial and biaxial) gives exactly σ/f<sub>y</sub>; stocky plate b/t = 20
  ≈ σ/f<sub>y</sub> (+1.4 %); pure shear √3τ/f<sub>y</sub>; b/t = 57 plate 0.415 vs RP-C201 0.420 / Faulkner 0.418.
* CSR proportion coefficients (100/125, 75/41/22, 12, b<sub>f</sub> ≥ 0.25 h<sub>w</sub>, √(235/R<sub>eH</sub>)) match IACS
  CSR BC&OT Pt 1 Ch 8 Sec 2.

### 6.2 Findings

**U1 — Ultimate capacity can exceed the global (column) elastic buckling load.** High · Theory
`solver.py:4913–4925` sets the Perry lateral-deformation moment to zero as soon as the load factor
passes the global elastic factor; `5015–5020` keeps the Ritz global curvature stresses
"diagnostic only"; `5781–5790` uses the elastic cut-off only when no yield is found. Load steps
can jump over λ<sub>E</sub>, Newton then converges on an unstable branch (reviewer: lowest tangent
eigenvalue −44†), and first yield is found near f<sub>y</sub>. Flat bar 120×10, L = 7 m, s = 800, t = 16,
σ<sub>x</sub> = 100: Euler stress with the full plate 20.5 MPa, so UF ≥ 4.88; the solver reports ultimate
UF 0.2885 (selected UF 0.33, valid, confidence "medium"). The *buckling* UF of the same run is
4.88, and RP-C201 gives 6.52. T 120×14/80×12 likewise: 0.29 vs ≥ 1.55. The reviewer found this in
24 of 432 scanned panels†. Since the API's default acceptance is `'ultimate'`, the selected UF passes.

**U2 — A "pressure-dominated" rule lifts the ultimate capacity to the elastic eigenvalue.** High · ECOSYSTEM_GUIDE
`solver.py:5809–5839`: when pressure alone produces ≥ 5 % yield utilisation
(`pressure_dominated_yield_preload_ratio = 0.05`), first-yield is dropped from the buckling-strength
envelope and `reported_ultimate = max(raw_ultimate, buckling_strength)` lifts the ultimate capacity to
the elastic factor. T 200×12/100×15, s = 500, t = 18, L = 1500, σ<sub>x</sub> = 200: ultimate UF 0.620 (p = 0),
0.638 (p = 0.05), **0.252** (p = 0.10 … 0.40) — below σ<sub>x</sub>/f<sub>y</sub> = 0.563. The reviewer suspects it exists
to satisfy the verification harness gate `max_uf_order_violations: 0` (ultimate UF ≤ buckling UF).

**U3 — Load factors and pressure side are ignored by the SemiAnalytical path.** High · Internal
`anystructure_panel_input` (`solver.py:927–998`) never reads `_stress_load_factor`,
`_lat_load_factor` or `_overpressure_side`; the RP-C201 path applies them. README panel with
`load_factor_stresses = load_factor_pressure = 1.5`: SemiAnalytical UF stays 1.0886, RP-C201 plate UF
0.282 → 0.423. Pressure is always taken as positive.

**U4 — Sniped stiffeners: ultimate UF decreases as lateral pressure increases.** High · Theory
`solver.py:5008–5014`: `pi_moment_stress = (-common_moment + M_snipe_PI)/Z_p` with
`common_moment = M_pressure + M_Perry`, so pressure and the Perry moment always *relieve* the
plate-induced branch. For a simply supported stiffener with pressure on the plate side the plate
is compressed by psL²/8 at mid-span, and both imperfection directions must be checked.
T 200×12/250×20 sniped, L = 3500, s = 800, t = 10, σ<sub>x</sub> = 100: UF 0.574 → 0.522 → 0.478 for
p = 0, 0.05, 0.10 MPa (with the sign corrected†: 0.612 → 0.678 → 0.755). Continuous stiffeners were
not affected in the reviewer's scan† because the stiffener-induced branch governs.

**U5 — Mixed-sign axial stress keeps only the compression end.** High · Internal
`_anystructure_axial_design_stress` (`solver.py:921–924`) returns `max(σx1, σx2)` when the signs
differ. σ<sub>x1</sub> = 20, σ<sub>x2</sub> = −300: ultimate UF 0.115 (identical to σ<sub>x1</sub> = σ<sub>x2</sub> = 20), although the
300 MPa tension end alone needs ≥ 0.845. The U3 path handles this correctly†. (The RP-C201 path has
the same `max` but also a von Mises check that catches the tension end.)

**U6 — CSR proportion check reads every dimension below 10 as metres.** High for the CSR flag (UF unaffected) · Internal
`_csr_dimension_mm` (`solver.py:890–896`) multiplies values < 10 by 1000; `FlatStru` passes mm.
t = t<sub>w</sub> = t<sub>f</sub> = 8 mm: net plate thickness used 7 998 mm, CSR vector `[1, 1, 1, 1]` ("green");
by hand s/t<sub>net</sub> = 700/6 = 116.7 > 81.4 → fail.

**U7 — Torsional deformation stress is dimensionally inconsistent.** Medium (unconservative) · Theory (source)
`solver.py:4660–4666`: E · y<sub>w</sub> · Φ<sub>0</sub> · (π/ℓ)² · (1/(1 − σ/σ<sub>ET</sub>) − 1) has units MPa/mm; the
CSR/DNV-CG-0128 form multiplies by the lever arm (t<sub>f</sub>/2 + h<sub>w</sub>). Reviewer: 0.18 MPa instead of
54.7 MPa at σ/σ<sub>ET</sub> = 0.9; restoring the lever arm raises UF by 2–7 %†. In effect the tripping
capacity becomes the elastic eigenvalue rounded to the load-step grid.

**U8 — Unexplained knock-down on the Perry moment.** Medium · ECOSYSTEM_GUIDE (source)
`solver.py:4890–4900`: factor 1 − λ⁴/12 for λ ≤ 1.56, else 3/(4λ) — discontinuous at 1.56
(0.507 vs 0.481) and without a cited source. Without it the solver matches an independent Perry
first-yield calculation (w<sub>0</sub> = L/1000) within 0.3 %; with it UF drops by up to 10 %†. It also
contributes to U1.

**U9 — Stress gradients in the yield check vary along the wrong direction.** Medium · Internal
S3 (`4731–4739`) varies σ<sub>y</sub> across the panel width (in y) and U3 (`4772–4788`) varies σ<sub>x</sub> along x —
i.e. each stress varies along its own direction, which violates in-plane equilibrium and
contradicts the package's own RP-C201 path (σ<sub>y</sub> varies along the stiffener length). Effect
+3–7 %† when swapped. Combined with S2 (elastic analysis uses the mean), pure in-plane bending is
rejected as `zero-variable-load` and the API returns `inf`.

**U10 — Domain gates and confidence label.** Medium · Source
Gates (`solver.py:247–252`) are wider than the PULS manual limits recorded in the same file
(`78–88`): plate b/t 250 vs 200; web 180 vs 90 (35 for flat bars); b<sub>f</sub>/t<sub>f</sub> 45 vs 15. Confidence is
"medium" for every valid run unless mode-convergence checking is enabled (`5307–5309`); the
conditions in U1/U2 and manual-limit violations never lower it.

**U11 — Low.** Panels that do not yield before λ = 100 take the elastic cut-off as capacity
(`5781–5790`), under-reporting UF for very light loads; linear interpolation across singular
jumps quantises capacity to the load step (`5124–5135`); silent fallbacks
(`except Exception: pass` at `4388`, `lstsq` on `LinAlgError` at `4337`, swallowed exceptions in the
batch predictor at `1266`) — ECOSYSTEM_GUIDE requires fallbacks to be explicit and reported; Perry
w<sub>0</sub> = L/1000 hard-coded (`4922, 4935`) ignoring the config value; CSR flange outstand measured
from the web face for T but from the centreline for angles (`1666–1669`); the CSR result metadata
embeds a local Windows file path as `source`.

**On calibration.** The constants listed in §5.3 and U2/U8 have no derivation in the code. The CLI
verification harness gates changes on error against PULS reference data
(`baseline_buckling_mae_delta`, `max_uf_order_violations`, `solver.py:6772–6890, 8600–8630`). This
repository's history is a single import commit, so whether the constants were tuned against PULS
cannot be established here; it should be documented either way, as the ECOSYSTEM_GUIDE requires.

## 7. Recommendations (from the first pass; see §0 for what has been done)

1. **Until U1/U2 are fixed, do not use the SemiAnalytical *ultimate* UF for acceptance.** (Done: U1, U2 fixed.) In every
   reproduced U1 case the SemiAnalytical *buckling* UF was the safe one (4.88 vs 0.29). Consider
   making `'buckling'` the default acceptance for that method, or capping the reported ultimate
   capacity at the global elastic factor, and flagging such runs as low confidence.
2. **Fix the internal defects first** — they are certain, local and do not depend on the wording of
   the standards: P1, P2, P5, P6, P7, P11, C1, C2, C3, C7, C8, C10–C14, U3, U5, U6, S4, S7.
3. **Confirm against your copy of the standards, then fix:** P3 (single z\*), P4 (I<sub>t</sub> factor),
   P8 (τ<sub>crg</sub> default), P9 (C<sub>τ</sub>), P10 (t<sub>w,red</sub>), C5 (column check for unstiffened shells),
   C6 (0.9/λ̄²), C9 (γ<sub>M</sub> scaling).
4. **Semi-analytical modelling:** choose the number of longitudinal half-waves from the aspect ratio
   (at least ⌈a/b⌉ + 2) instead of a fixed 3 (S1); build the geometric stiffness from the actual
   linear stress distribution (S2, U9); keep the column (Perry) term beyond λ<sub>E</sub> and check
   tangent-stiffness positive-definiteness, or cap capacity at λ<sub>E</sub> (U1); remove the pressure lift
   (U2); correct the plate-induced branch sign (U4) and the torsional lever arm (U7); either derive and
   document or remove the empirical constants (S3, S5, S6, S8, U8).
5. **Tests.** Replace parity-only golden values with independent references: closed-form plate
   buckling coefficients, Euler column, hand-calculated RP-C201/RP-C202 cases; plus invariance and
   monotonicity properties that would have caught several findings directly — σ<sub>y1</sub> ↔ σ<sub>y2</sub> relabelling
   (P11), sign of σ<sub>m</sub> on a full cylinder (C4), UF non-decreasing in pressure and in load factor
   (U2, U3, U4), ultimate capacity ≤ global elastic capacity for column-type failure (U1), and
   continuity across formula branch limits (P9, P10, C3, C6, C9). Every defect fix should land with a
   regression test that fails without it (ECOSYSTEM_GUIDE §3). Correct
   `tests/test_structure.py::test_plastic_sec_mod` together with P13.
