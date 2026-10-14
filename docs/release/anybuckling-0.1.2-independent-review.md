# ANYbuckling 0.1.2 — independent release review

- **Reviewer:** independent Claude subagent (Claude Code, model claude-opus-5-5), spawned in session
  https://claude.ai/code/session_01KGRyBUB2uYou3Xyq1TXwMM. The reviewer did not build the artifacts or write the
  evidence file. The repository was not modified. Nothing was committed, pushed, tagged or installed into the
  shared environment.
- **Date:** 2026-10-14
- **Environment:** Linux 6.18.44 x86_64 (Claude Code cloud container), Python 3.13.16, fresh venv
  `review/venv` (pip 26.2.1): wheel + numpy 2.5.4, scipy 1.18.1, pytest 9.1.1 from PyPI, numba absent.
- **Scratch directory:** `/tmp/claude-0/-home-user-ANYbuckling/8ccb8bb8-454e-5a0c-96e0-0d63995c8e08/scratchpad/review`
  (`src/` = `git archive f123816`, `old/` = `git archive v0.1.1`, `whl/`, `sdist/` = extracted artifacts,
  `scripts/` = review scripts).

## Inputs

| Input | Identity |
|---|---|
| Artifact-source commit | `f123816703ccb779b24b57bc22dd881229464822`, tree `b88b0cd89139edd856531b38696b4e37453e3198` (both confirmed with `git rev-parse`) |
| Previous release | tag `v0.1.1` → `a871d5a3c466666b79f3ce3a015a2cfd7534376b` |
| `anybuckling-0.1.2-py3-none-any.whl` | 165072 bytes, sha256 `b22a10ee79e773f479f19788f006daca0d7e188ad3c64e9309122749b00a05a5` |
| `anybuckling-0.1.2.tar.gz` | 195663 bytes, sha256 `8f68e0e5a549b6774c12e8e1bf0bd764f68e303c52ae43e2dc877986e07f38a7` |
| `SHA256SUMS` | 192 bytes, sha256 `55eb9636da5d58c7ac332a814bd340f0acb7b6f2135355365970e8fe53dccf75` |
| Qualification evidence `anybuckling-0.1.2-qualification.json` | sha256 `9b106c711e9427352c4c486d80456e765eba8688caacde247ed22e4cb5ef4df4` |
| `tools/verify_release_authority.py` @f123816 | sha256 `25aec5566bb0e4ebfbd7c1005c35bcd180dc58f2478008517bcae01643596033` |
| `.github/workflows/publish-release-assets.yml` @f123816 | sha256 `1471d3dae4b3ab4fbcce03cdd46b194011dd3e14c3edb1883e44509f12bbc970` |
| `docs/release/anybuckling-0.1.2-notes.md` @f123816 | sha256 `853c41b48169c31508521d63ddbdf4528a2fcdee8bb569735da1b4acbc3a814c` |
| `ECOSYSTEM_GUIDE.md` @f123816 | sha256 `a9827dfb710c2bbcedbb6516105438b8defa0cb6b8dd2279cf26812ff921ebdc` |

## Check 1 — Artifact identity: PASS

Commands: `sha256sum *.whl *.tar.gz`, `stat -c '%n %s'`, `od -c SHA256SUMS`, then a Python check that
`SHA256SUMS` equals `"".join(f"{sha256.lower()}  {name}\n" for name in sorted(names))`, the exact form built by
`verify()` in `tools/verify_release_authority.py`.

- Sizes and digests match the evidence file and `SHA256SUMS` (values in the table above).
- `SHA256SUMS`: lowercase hex, two spaces, names sorted (`.whl` before `.tar.gz`), LF endings, trailing newline.
  It is byte-equal to the verifier's expected form: **True**.
- The verifier's own `_verify_wheel` and `_verify_sdist`, imported from the f123816 copy, pass on both files.
  They check for no duplicate, unsafe or link members, exactly one METADATA/PKG-INFO, and Name/Version
  ANYbuckling/0.1.2.
- Note for the ledger step: the verifier's `SHA256_RE` is `[0-9A-F]{64}` and it compares
  `hexdigest().upper()`, so the ledger `sha256` values must be **uppercase**, as they are in the 0.1.1 ledger.
  The evidence file uses lowercase. That is fine for evidence, but the values must not be copied verbatim into
  the ledger.

## Check 2 — Artifact content: PASS

Commands: `git archive f123816 | tar -x -C src`, `unzip` the wheel, `tar -xzf` the sdist, then `cmp` each file
against the same path in `src/`.

- Wheel: all 15 files under `anybuckling/` are byte-identical to `src/src/anybuckling/*`. The file set equals the
  package file set in git (`diff` of sorted lists: identical). `dist-info/licenses/LICENSE` is identical to
  `LICENSE`. Every RECORD hash and size was recomputed: 0 mismatches.
- Sdist: all 32 non-generated files (`src/anybuckling/**`, `tests/*.py`, `README.md`, `LICENSE`,
  `pyproject.toml`) are byte-identical to git. The only other files are generated: `PKG-INFO`,
  `src/ANYbuckling.egg-info/*`, and `setup.cfg`, which contains only the standard setuptools
  `[egg_info] tag_build = / tag_date = 0` stanza. The members are 9 directories and 39 regular files, with no
  links.
- Metadata: wheel METADATA and sdist PKG-INFO both have `Name: ANYbuckling`, `Version: 0.1.2`,
  `Requires-Python: >=3.10`, `Requires-Dist: numpy, scipy`, extras `excel` (xlwings) and `dev`.
- No unexpected files: no PDFs, no review documents (`docs/` is not shipped), no caches or `__pycache__`, no
  `.git`. A grep for key/secret/password markers matched only the GPL text in `LICENSE`.

## Check 3 — Version agreement at f123816: PASS

| File | Value |
|---|---|
| `pyproject.toml` | `version = "0.1.2"` |
| `src/anybuckling/__init__.py` | `__version__ = "0.1.2"` (also confirmed at runtime from the installed wheel) |
| `.github/workflows/publish-release-assets.yml` | `--ledger docs/release/anybuckling-0.1.2-ledger.json`, `--expected-terminal ACCEPTED_ANYBUCKLING_0_1_2_RELEASE`, `--distribution ANYbuckling`, `--version 0.1.2`, `--artifact anybuckling-0.1.2-py3-none-any.whl`, `--artifact anybuckling-0.1.2.tar.gz` |
| `tests/test_release_workflow.py` | `VERSION = "0.1.2"`, `EXPECTED_TERMINAL = "ACCEPTED_ANYBUCKLING_0_1_2_RELEASE"`, `WRONG_TAG = "v0.1.1"` |
| `docs/release/anybuckling-0.1.2-notes.md` | "ANYbuckling 0.1.2 — release notes" |
| Artifacts | `0.1.2` in file names and metadata |

The bump 0.1.1 → 0.1.2 is patch-only (guide section 5). The tag to be created is `v0.1.2`. The verifier enforces
`tag == f"v{version}"`.

## Check 4 — Tests against the installed wheel: PASS

```
python3 -m venv review/venv
review/venv/bin/pip install rel/dist/anybuckling-0.1.2-py3-none-any.whl pytest
cd review/src && ../venv/bin/python -m pytest -q -p no:cacheprovider tests
→ 252 passed in 10.06s
```

- The tests come from `git archive f123816`, and the package is installed from the wheel.
- Import location: no `sys.path` manipulation exists in `tests/` or `pyproject.toml`. A probe test run in the
  same pytest session printed
  `ORIGIN .../review/venv/lib/python3.13/site-packages/anybuckling/__init__.py`, so the package was imported from
  site-packages, not from `src/`.
- This reproduces the evidence claim "252 passed" (wheel, numba absent). I did not re-run the numba or sdist
  variants.

## Check 5 — Public API backward compatibility v0.1.1 → f123816: PASS (no incompatible change)

Method: `scripts/apicmp.py` does an AST comparison over every module in `src/anybuckling`. It covers public
top-level functions, classes, constants and `__init__` re-exports, each class's public methods and `__init__`,
parameter names, positions, required/optional status and kw-only status, and dataclass fields with or without
defaults. `scripts/defaults.py` compares default values.

- 603 public items in v0.1.1 and 621 at f123816. **0 removed names, 0 removed or renamed parameters,
  0 parameters moved or made required, 0 new required parameters or fields.**
- Added items, all additive:
  - `FlatStru.set_girder(span=None)`
  - `helper_cylinder_stress_to_force_to_stress(length_unit=None)`
  - `CylinderAndCurvedPlate.unstiffened_column_buckling`
  - `S3SolverConfig.auto_mode_count` (with default)
  - three `MembraneField` fields (with defaults)
  - new solver helpers and constants (`dnv_rp_c201_plate_buckling_coefficient`, `axial_buckling_driver`,
    `transverse_buckling_driver`, `stiffener_torsion_constant`, `plate_mode_tuples`,
    `applied_plate_stress_field`, `anystructure_panel_inputs`, `anystructure_adapter_diagnostics`,
    `CLAMPED_/FREE_EDGE_ROTATIONAL_SUPPORT`, `GLOBAL_*`)
- Changed defaults:
  - `CylStru.set_end_cap_pressure_included_in_stress(is_included)` True → False. Runtime check
    (`scripts/spot.py`): in 0.1.1, `default == True == False` gave identical results (ψ = 2 was always used). In
    0.1.2, `default == False` gives the same ψ = 2 path, and `True` now gives ψ = 4 and warns. The default call
    is therefore unchanged. Callers who explicitly passed `True` now get different results, with a
    `UserWarning`. This is documented in the notes as a defect fix.
  - `S3SolverConfig.include_pressure_dominated_yield_in_buckling_strength` False → True. This is the U2′ fix:
    the notes say "the buckling UF is never below the ultimate UF". The notes do not name the flag, but the
    behaviour is documented and the flag is listed in `docs/review/CALCULATION_REVIEW.md`.
- Result dictionary keys (runtime, `scripts/keys.py`, same inputs on 0.1.1 source and the 0.1.2 wheel):
  - 3 flat domains × 2 methods and 8 cylinder domains were compared.
  - No key was removed in any case.
  - One key was added for "Unstiffened shell" (the new column check).
- Unit-semantics changes, signature-compatible and documented as fixes:
  - `set_length_between_girder` now converts mm → m.
  - `helper_cylinder_stress_to_force_to_stress` infers the unit when `length_unit=None`.
  - `set_girder` without `span` keeps the 10 m placeholder, warns, and gives results identical to
    `span=10000` (verified at runtime).
- Additional evidence: the v0.1.1 test suite was run against the 0.1.2 wheel. Result: 98 passed and 6 failed.
  All 6 failures are pinned numerical values. Each has an updated value and justification in the 0.1.2 tests,
  matching documented fixes:
  - Table 3-3 torsion 5.34: 2 cylinder tests
  - P13 plastic section modulus: 1
  - S1/S5 semi-analytical: 3
  None of the failures is an API error.

## Check 6 — Release-notes accuracy (spot checks): PASS

`git log --oneline v0.1.1..f123816` lists 23 commits: plate, cylinder and semi-analytical fixes, review docs,
S5, performance, and the release prep. Claims checked:

1. **"`FlatStru.set_girder(..., span=L_G_mm)` sets it. Omitting `span` keeps the previous 10 m placeholder,
   with a warning (P1)."** Confirmed:
   - The code is in the `api.py` diff (commit 23068b9).
   - Runtime: no span gives a `UserWarning`; `span=10000` gives no warning; the results of the two are
     identical.
2. **"V ≥ V_Rd now fails instead of raising (P10)."** Confirmed at runtime with the P10 regression case (p = 0.20,
   FB 200×6):
   - 0.1.1 raises `ZeroDivisionError`.
   - 0.1.2 returns Shear capacity UF 1.049 (> 1).
   - Minor wording point: the notes say "UF = ∞", but the reported shear UF here is finite (> 1). Failure is
     reported either way.
3. **"N_Sd ≥ N_E now fails (UF = ∞) (P17)."** Commit 793fc5b adds `math.inf` returns in `plates.py` and a
   regression test. That test passes (Check 4).
4. **"`set_end_cap_pressure_included_in_stress(is_included=False)` now does what its name says. The default
   keeps the previous results; `True` (ψ = 4) warns."** Confirmed at runtime (see Check 5).
5. **"`set_length_between_girder` stores metres."** Confirmed in the diff: `api_helpers.mm_to_m(val)`.
6. **"Global stiffness: H = D + G I_t/(2s) instead of Huber (S5)."** Confirmed:
   - Commit 26cbbc4.
   - `solver.py` docstring and code around `stiffener_torsion_constant` (lines ~2082–2107).
7. **"Supports: free edges are rejected (S4)."** Confirmed: `FREE_EDGE_ROTATIONAL_SUPPORT = "FS"`, rejected at
   `solver.py` ~836–840 and ~2706.
8. **"Column buckling: unstiffened shells now get the column check."** Confirmed:
   - New method `unstiffened_column_buckling`.
   - New result key for "Unstiffened shell".
9. **"`S3SolverConfig(auto_mode_count=False)` restores the old basis."** The field exists with default True. I
   did not verify numerically that it reproduces 0.1.1.
10. **"Table 3-3 torsion 5.34"** and **"Plastic section modulus corrected (P13)"**: confirmed by the changed
    pinned values in `tests/test_cylinders.py` and `tests/test_structure.py`. The old values fail against
    0.1.2, as expected.
11. **Performance, "about 28 ms (19 ms in 0.1.1)":** not reproduced on my panel. The reference stiffened panel
    took 25.9 ms per call on 0.1.1 and 25.6 ms per call on 0.1.2 (numba absent). The notes' panel is not
    specified, so this is informational only. Performance evidence is not correctness evidence (guide
    section 4).

## Findings

| # | Severity | Finding |
|---|---|---|
| F1 | Info (process) | f123816 is not yet on `origin/main`, which is at its parent `7084089`. The verifier requires the ledger commit to be reachable from `refs/remotes/origin/main` and to be the sole child of f123816, so f123816 must be on main before the ledger and tag are created. |
| F2 | Info (process) | The ledger `sha256` fields must be uppercase hex (`SHA256_RE = [0-9A-F]{64}`). The evidence file's lowercase digests must be uppercased when the ledger is written. |
| F3 | Low (docs) | The notes say V ≥ V_Rd gives "UF = ∞". In the P10 regression case the shear UF is finite (1.049 > 1). The failure is still reported, which is the substantive claim. |
| F4 | Low (docs) | The `S3SolverConfig.include_pressure_dominated_yield_in_buckling_strength` default flip (False → True) is described in the notes only by its effect (U2), not by name. |
| F5 | Info | The performance figure (28 ms vs 19 ms) was not reproduced; my reference panel ran at ~26 ms on both versions. Not a release gate. |
| F6 | Info (guide section 3, outside the requested scope) | The repository has no CI test workflow; the evidence file records this as a limitation. The PULS benchmark gates were not run. Both are recorded in the evidence `limitations`. |

None of these findings is a failed check or a backward-incompatible public API change.

VERDICT: ACCEPT
