# JMP feature parity work

The official JMP/JMP Pro 19.2 feature inventory and this bounded implementation
batch are recorded. **Full JMP parity remains unfinished.** The authoritative
backlog is [jmp-feature-matrix.json](jmp-feature-matrix.json), with scope and
primary sources in [JMP_PARITY.md](../JMP_PARITY.md).

The inventory contains 552 heterogeneous feature/method/option rows: 287 partial,
197 absent, 37 provided by the host database, 6 earlier verified limited contracts,
17 GUI exclusions, 0 output/scope questions, 4 external plugins and 4 separate
product entries. These counts do not constitute a completion percentage. A
20.0 help URL returns the same index as 19.2; a distinct 20.0 feature delta has not
been verified. No installed JMP was used for numerical comparisons.

## Verified multivariate PLS checkpoint (2026-10-24)

PLS2 and indicator PLS-DA fit/predict add four operations, for **277 analytics
functions /428 overloads**. Both NIPALS/SVD and SIMPLS include VIP; regression
prediction includes explicit normality-based conditional-score approximate
mean/individual Student limits. Four feature rows move from absent to partial.
The matching macOS arm64 dev/Debug source passes1,391 Rust tests in69 suites,
0 failures,2 existing ignored;20 SQL cases/626 assertions,0 skips;428 catalog
examples;strict Clippy and formatting. Independent NumPy/SciPy/sklearn checks
pass7,220 numerical fields and34 malformed NULL requests across C ABI and SQL.
[Contracts](jmp-pls-api.json), [validation](evidence/jmp-pls-validation.json).
Full JMP completion remains open.

## Verified zero-inflation checkpoint (2026-10-24)

Five natural-support unpenalized shared-inflation regression fits and portable
prediction add six operations: **273 analytics functions /424 total overloads**.
The new feature row is partial; whole-family and full JMP completion remain open.
The matching dev/Debug CPU build passes1,381 Rust tests in68 suites,0 failures,
2 existing ignored;19 SQL cases/616 assertions,0 skips;424 catalog examples;
strict all-target analytics Clippy and scoped formatting. Independent SciPy and
mpmath references check coefficients, full joint covariance/HC0, saved prediction
and numerical tails at the actual C ABI and loaded SQL boundaries.
See [contracts](jmp-zero-inflation-api.json),
[validation](evidence/jmp-zero-inflation-validation.json) and the
[active continuation record](jmp-completion-progress-20261023.json).

## Existing API mapping audit (2026-10-23)

Four rows previously marked absent now link to existing bounded implementations:
`glmm_binomial_response`, `glmm_count_response`, `sem_fit_indices`, and
`text_analysis`. Twelve API entries gained explicit feature links. The source
hashes for GLMM, SEM and text match their recorded native SQL validation
checkpoints. Each corrected row remains partial and lists unsupported contracts.
See [the audit record](evidence/jmp-feature-mapping-audit-20261023.json).

At this historical audit, the registered API counts remained 267 analytics functions and 418 total overloads.
No algorithms or APIs were added, and no runtime tests were rerun for this
metadata/documentation correction. The previous runtime evidence and its original
status counts remain unchanged. This audit corrects four confirmed omissions;
it does not certify every other absent row as having no bounded implementation.

The following sections describe historical checkpoints, including their then
current counts and work status. They do not change the scope of this mapping audit.

## Historical until-achieved continuation (2026-10-22)

The user explicitly requested continuing until the objective is achieved. The
persistent goal remains active; the batches below are verification checkpoints,
not overall completion. GUI operations and GPU dependencies are excluded; CPU
numerical methods and their usable SQL interfaces remain in scope.

The previous verified checkpoint added 62 operations to the earlier 119, giving
181 analytics functions and 332 total extension overloads. It passed 1,099 Rust
tests in 48 suites (2 existing external-data tests ignored), all 9 SQL cases with
394 assertions, all 332 catalog examples, strict Clippy and scoped formatting.
The 62-function SQL checker validated 47 operations with numerical fixtures,
107 result fields, 11 saved-model round trips and 16 malformed NULL requests.
See [jmp-continuation-validation.json](evidence/jmp-continuation-validation.json).
The README now contains all 209 individual function examples and exact limits.

Current inventory: 235 partial, 205 absent, 37 host-provided, 6 earlier verified
limited contracts, 17 GUI exclusions, 44 output/scope entries, 4 external plugins,
4 separate products (552 heterogeneous rows; no completion percentage).

The latest checkpoint verifies the next 28 operations, giving 209 analytics
functions and 360 total overloads; its detailed results are below. New isolated
work is active on IRT, SEM/CFA and degradation/stability/cumulative exposure.
These candidate methods are not counted as registered or verified. The single task-owned build
tree is retained for this immediate continuation; the prior batch tree had been
removed after its evidence was saved.

## Earlier delivered implementation batch

Added 95 documented CPU SQL functions to the original 24, for 119 analytics APIs
and 270 total extension overloads. The complete function-by-function input,
ordered parameters, SQL examples, model reuse, limitations and remaining gaps
are in the [README](../README.md#cpu-analytics-complete-sql-api-reference).

- Inference (13), OLS/GLM/regularization (6), quality/MSA (18).
- Resampling/simulation/planning/effect sizes (11), distributions/lifetime (12).
- Time series (10), multivariate/classification/clustering (8).
- ML/REML random-intercept fit/prediction (2), Cox fit/prediction (3), CSS ARIMA (1).
- CART/forest/first-order boosting/single-hidden-layer CPU neural models (4).
- Box-Behnken/LHS/mixture/fractional/Hadamard designs, randomization/evaluation (7).

No GUI or GPU work and no new runtime dependencies were introduced. Existing
user work in other worktrees was preserved. Ryer tools were not exposed; indexed
Tsubame discovery was used. Changes remain uncommitted on
`feat/sql-analytics-20261021`, based on `985af99`.

## Verification

Behavioral RED was observed before implementing new methods and numerical fixes,
followed by GREEN, independent analytical/NumPy/SciPy fixtures, invariants and
boundary cases. Public dispatch also received RED/GREEN integration coverage.

- Rust workspace: 1,006 passed, 0 failed, 2 pre-existing external-data tests ignored
  across 40 suites. Scoped formatting and strict analytics Clippy passed.
- Native macOS arm64 DuckDB 1.5.3 dynamically loaded the extension. C++ Debug
  AddressSanitizer/UndefinedBehaviorSanitizer instrumentation was enabled; Rust
  libraries were not sanitizer-instrumented. All 8 SQL test cases passed.
- All 270 catalog examples executed with non-NULL results and complete metadata.
- Original SQL checker: 25 numerical/model checks, 7 invalid inputs, selected
  vectors and actual molecular descriptor integration passed.
- New SQL checker: all 95 additional operations, 199 numerical fields, 3 saved
  model SQL round trips, 14 invalid inputs, 7 models ignoring suffix labels,
  7 numerical regression cases and 12 selected/dynamic-vector rows passed.
- Separate reference evidence records model/inference/time-series/Cox,
  distribution/ARIMA and design comparisons with reproducible harnesses.

Commands, source/artifact hashes, environment and limits are preserved in
[jmp-sql-validation.json](evidence/jmp-sql-validation.json). The single task-owned
temporary build tree is removed after final documentation verification.

## Remaining work

Implement and verify the absent/partial numerical contracts in the inventory,
including GLMM, random slopes/repeated covariance, SEM, functional/spectral
analysis, advanced DOE, nonlinear and consumer/text/genomic methods. Existing
family-level entries remain partial even when individual methods pass. Windows,
Linux and Wasm integration for this expansion, and actual JMP numerical
concordance, remain unverified.


## Verified 209-function checkpoint (2026-10-22)

Added28 operations across mixed_advanced,glmm,causal,omics,design_classical,market. Current209analytics/360total overloads;1167Rust tests in54suites GREEN,2existingignored;11nativeSQLcases/420assertions GREEN;all360catalogexamples nonNULL. See [verification](evidence/jmp-models-validation.json) and complete README contracts. Current552-row inventory: partial 227, database_existing 37, absent 213, external_plugin 4, clarification_pending 44, implemented_verified 6, excluded_gui 17, separate_product 4. The extra row distinguishes classical item reliability from IRT Item Analysis. Earlier summaries above describe their historical snapshots.

The until-achieved goal remains ACTIVE. New isolated TDD work proceeds on IRT,SEM/CFA and degradation/stability/cumulative exposure;these are not counted as registered/verified in this checkpoint. All other absent/partial contracts in the matrix remain required work.

The next isolated stage also includes supplied-kinship Gaussian ML/REML and saved cross-kernel BLUP. Its 11 candidate tests cover exact balanced component estimates, original-design REML units, saved predictions, flat identification, near-endpoint optima, perfect-fit and physical-underflow rejection. It is outside the209 compiled snapshot until its public/native gates pass.


## Verified225-function checkpoint (2026-10-22)

Added16 operations across kinship,IRT,SEM and degradation. 1236Rust tests/58suites pass,2existingignored; nativeSQL 12cases/430assertions pass;all376catalogexamples executed. [Current evidence](evidence/jmp-specialized-validation.json). Current552-row status: partial 235, database_existing 37, absent 205, external_plugin 4, clarification_pending 44, implemented_verified 6, excluded_gui 17, separate_product 4. The until-achieved goal remains ACTIVE; new time-series/state-space and further functional preprocessing contracts are the next isolated stages. Earlier paragraphs describe their historical checkpoints.


## Final254-function checkpoint; push and stop (2026-10-22)

Current29operations completed and matching native verified:1302Rust tests/63suites,13SQLcases/439assertions,405catalogexamples,29publicC ABI fixtures,9SQL savedmodel/reference roundtrips. [Final evidence](evidence/jmp-streams-validation.json). Current552-row inventory: partial 254, database_existing 37, absent 186, external_plugin 4, clarification_pending 44, implemented_verified 6, excluded_gui 17, separate_product 4. These are bounded partial contracts;overall JMP parity remains incomplete. The user explicitly requested push then stop,so no further feature stage is included. Previous ACTIVE and continue-until-achieved entries above describe the earlier instruction;the current stop instruction supersedes further autonomous implementation after this checkpoint is pushed.
