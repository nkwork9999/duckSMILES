# SQL analytics implementation

This record covers the original 24-function revision. The subsequent 95-function
expansion and its current verification are recorded in
[jmp-parity-progress.md](jmp-parity-progress.md); the complete 119-function
reference is in [README.md](../README.md#cpu-analytics-complete-sql-api-reference).

Worktree: `ducksmiles-analytics-20261021`, branch `feat/sql-analytics-20261021`.
Starting revision: `985af99` (origin/main). Existing worktrees left untouched.

Authorized scope: non-GUI analytics discussed with the user, integrated in order:

1. PCA and deterministic k-means.
2. Ridge/PLS1 regression, prediction and fold-local cross-validation.
3. Scaffold/group-controlled fold assignment.
4. Factorial and central-composite DOE, candidate-constrained D-optimal selection,
   and quadratic response surfaces.
5. Gaussian-process candidate scoring, expected improvement, uncertainty and
   diversity selection for active learning / Bayesian optimization.

Additional non-GUI categories from the comparison: descriptive statistics,
inferential tests and multiplicity, process capability/control limits, and
explicitly scoped mixed-model / survival analysis.

Architecture: a separate dependency-free Rust analytics crate linked into the
existing DuckDB extension. SQL receives typed numeric lists and returns
versioned JSON. Reusable models expose parameters and preprocessing metadata.
Null/malformed/nonfinite/oversized inputs must not panic or cross FFI unwinding.
All learning preprocessing must be fitted inside each training fold.

Verification: tests first with independently derived fixtures, Rust workspace
suite, focused lint/format checks, DuckDB integration and function-catalog
examples. Record exact supported methods and limits rather than claiming JMP
feature equivalence. No GUI changes or publishing requested.

## Progress

- Repository identified; no Ryer tools exposed in this session.
- Dedicated worktree created; baseline workspace tests started.
- All five stages integrated, plus the explicitly scoped statistical methods.
- Added 24 documented SQL functions and 54 Rust regression tests.
- Final Rust workspace: 865 passed, 0 failed, 2 existing ignored tests.
- Native SQL extension suite: 7 cases, 367 assertions, all passed, no skips.
- Function catalog: all 175 overloads documented; all 175 examples executed.
- SQL model/boundary checker: 25 numerical/model checks, 7 malformed inputs,
  dictionary/selected-vector rows and molecular-descriptor integration passed.
- NumPy/SciPy reference checks: 18 seeded fixtures; after correcting Gaussian
  tail cancellation, GP means/variances/EI agree within 1.12e-15, including exact
  candidate ranking. Other algorithms agree within approximately 1.3e-14.
- All three documentation SQL blocks executed successfully, including saved
  linear/quadratic model prediction and the end-to-end molecular workflow.
- Strict analytics Clippy, scoped Rust/C++ formatting and diff checks passed.
- Validation evidence: `tasks/evidence/sql-analytics-validation.json` and
  `tasks/evidence/sql-analytics-reference.json`.
- Temporary build products removed after verification; implementation, tests,
  examples and evidence retained on the dedicated local branch.

## Scope limits

This delivers the requested calculation workflow, not full JMP parity. D-optimal
design is candidate-constrained local exchange. Bayesian optimization ranks a
finite candidate pool with explicit hyperparameters. Mixed effects are a balanced
one-way random-intercept ANOVA estimator; survival analysis is Kaplan-Meier.
PCA rejects extreme offsets that prevent representable centering. The native SQL
boundary was verified on macOS arm64; other platform builds remain unverified.
