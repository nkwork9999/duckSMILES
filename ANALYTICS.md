# SQL analytics

This guide describes the original 24-function foundation. The complete current
279-function API, individual examples, numerical limits and validation results
are in [README.md](README.md#cpu-analytics-complete-sql-api-reference).
The researched JMP feature backlog is in [JMP_PARITY.md](JMP_PARITY.md).

The analytics API adds numerical exploration, regression, molecule-aware evaluation,
experimental design, candidate selection, and a defined set of statistical methods
to duckSMILES. Computation runs in the extension's dependency-free Rust core. It
has no GUI and does not claim the full statistical coverage of JMP/JMP Pro.

These functions require a build containing the analytics module; an older
community extension package may not contain them. Check the loaded build:

```sql
SELECT function_name, parameters, description
FROM duckdb_functions()
WHERE starts_with(function_name, 'chem_')
ORDER BY function_name;
```

## Common contract

All numeric functions have the SQL signature
`chem_function(data DOUBLE[][], parameters DOUBLE[]) -> VARCHAR`. They return JSON
with `"version":1`. The generic entry point
`chem_analyze(operation VARCHAR, data DOUBLE[][], parameters DOUBLE[])` calls the
same operations, using the operation names in the tables below. The separate
`chem_scaffold_folds(smiles VARCHAR[], parameters DOUBLE[])` wrapper computes
molecular scaffolds before group partitioning; it is not a numeric
`chem_analyze` operation.

Each inner list is one observation. Its columns have the positional meaning
specified below. Use `[]::DOUBLE[]` for an empty parameter list. Optional numeric
parameters can be omitted from the end of a list. All input values must be finite,
all rows must have the same nonzero width, and lists cannot contain `NULL`.
Malformed, undefined, numerically singular, overflowing, or resource-limited
requests return SQL `NULL` for that request. Missing-value imputation and
categorical encoding are caller responsibilities. Individual JSON fields can be
`null` where explicitly documented; that differs from a failed SQL request.

JSON result indices and labels are zero based. DuckDB SQL list indexing is one
based. Returned rows, labels, and prediction arrays preserve input order unless a
result explicitly represents a ranking or a design. Always aggregate using a
stable unique key, for example `list([mw, logp] ORDER BY id)`. Save the feature
names, order, units, model JSON, input version, seed, and hyperparameters together
when retaining a model.

## Exploration

| SQL function / operation | Data rows | Parameters in order | Main JSON fields |
| --- | --- | --- | --- |
| `chem_pca` / `pca` | `[feature_0, ...]` | `[components=min(p,n-1), standardize=0]` | `mean`, `scale`, `components`, `scores`, `explained_variance`, `explained_variance_ratio` |
| `chem_kmeans` / `kmeans` | `[feature_0, ...]` | `[clusters=min(2,n), max_iterations=300, seed=0]` | `labels`, `centers`, `inertia`, `iterations`, `converged` |

PCA requires at least two observations and `1 <= components <= min(p,n-1)`.
It centers every feature and diagonalizes its sample covariance, using denominator
`n-1`. Set `standardize=1` to divide features by their sample standard deviations;
`0` uses their original units. Constant features use scale `1`. Component vectors
are rows of `components`, ordered by descending variance. The largest absolute
loading is assigned a positive sign for deterministic output. Explained-variance
ratios are zero for a wholly constant matrix. Repeated eigenvalues can have
nonunique component directions, so compare their subspaces when appropriate.

For a new row `x`, reuse the saved model without refitting:

```text
z[j] = (x[j] - mean[j]) / scale[j]
score[k] = sum_j z[j] * components[k][j]
reconstruction[j] = mean[j] + scale[j] * sum_k score[k] * components[k][j]
```

The reconstruction is truncated if only some components were retained. There is
no dedicated PCA prediction function; these saved arrays specify the transform.
Do not fit PCA on held-out rows before evaluating a downstream predictor.

K-means uses seeded k-means++ initialization followed by Euclidean assignment and
mean updates. Features are used in their supplied units: standardize or weight
them explicitly when needed. Inertia is the sum of squared distances to the
returned centers. Empty clusters are repaired, including when repeated points
produce fewer distinct positions than requested clusters. A capped run returns a
consistent partition with `converged=false`; this is not a global-optimum claim.
Cluster count is at most `min(n,256)`, and iterations must be in `1..10000`.

## Regression and evaluation

| SQL function / operation | Data rows | Required parameters in order | Main JSON fields |
| --- | --- | --- | --- |
| `chem_ridge` / `ridge` | `[features..., target]` | `[lambda]`, nonnegative | `intercept`, `coefficients`, `stable_predict_parameters`, `feature_means`, `feature_scales`, `metrics` |
| `chem_pls` / `pls` | `[features..., target]` | `[components]`, integer `1..p` | Same model fields, plus `components_used` |
| `chem_predict` / `predict` | `[features...]` | Saved `stable_predict_parameters`, or `[intercept, coefficient_0, ...]` | `predictions` |
| `chem_cross_validate` / `cross_validate` | `[features..., target, fold_id]` | `[method, hyperparameter]`; method `0`=Ridge, `1`=PLS1 | `out_of_fold_predictions`, `rmse`, `r2`, `fold_metrics`, `preprocessing` |
| `chem_group_folds` / `group_folds` | `[group_id]` | `[n_folds, seed]` | `folds`, `fold_sizes`, `n_groups` |
| `chem_scaffold_folds` | List of SMILES strings | `[n_folds, seed]` | Same fold result |

Ridge and PLS1 require at least two training rows. Both center the target and
center/scale each feature using the training population standard deviation
`scale[j] = sqrt(sum_i (x[i,j]-mean[j])^2 / n)`; constant columns use scale `1`.
Ridge minimizes

```text
sum_i (y[i] - intercept - sum_j beta[j]*z[i,j])^2
    + lambda * sum_j beta[j]^2
```

The intercept is unpenalized; the penalty applies to coefficients in standardized
feature units, and is not divided by the number of observations. Ridge uses an
augmented QR solve. Zero-penalty unidentified designs are rejected; use a positive
penalty or an appropriate reduced model. PLS1 supports a single response, deflates
predictive components, and can use fewer components than requested when predictive
rank is exhausted. Constant responses can produce zero predictive components.

Returned linear `coefficients` and `intercept` are converted to the original
feature units. `chem_predict` accepts raw feature rows; do not apply the stored
standardization a second time. For portable models, pass the saved
`stable_predict_parameters` vector unchanged. It contains
`[response_origin, response_mean_offset, feature_origins..., feature_mean_offsets..., coefficients...]`
and has length `3p+2`. The equivalent evaluation is
`response_origin + (response_mean_offset + sum_j coefficients[j]*((x[j]-feature_origins[j])-feature_mean_offsets[j]))`.
This preserves small differences around large numeric offsets. The older
`[intercept, coefficients...]` form remains accepted but can lose precision through
cancellation at extreme offsets. For example:

```sql
WITH model AS (
    SELECT chem_ridge([[0.0,1.0],[1.0,3.0],[2.0,5.0]], [1.0]) AS j
)
SELECT chem_predict(
    [[3.0],[4.0]],
    json_extract(j, '$.stable_predict_parameters')::DOUBLE[]
)
FROM model;
```

Fit `metrics` are in-sample RMSE and R². The linear-regression/CV R² field is JSON
`null` for a constant evaluation target. Cross-validation instead fits one model
per held-out fold and reports predictions in original row order, pooled RMSE/R²,
and each fold's train/test counts and metrics. Centering and scaling are fitted
only on each training fold (`preprocessing="training_fold_only"`). Each training
partition still needs at least two rows and an identifiable chosen model.
Hyperparameters are supplied explicitly; this is not nested tuning, automatic
feature selection, or automatic preprocessing of arbitrary external features.

Fold and group identifiers are exact nonnegative integers. There must be
`2..1000` nonempty folds. `chem_group_folds` keeps each group indivisible, assigns
larger groups first, and uses the seed to break ties while balancing row counts.
The number of distinct groups must be at least the requested number of folds;
exact balance is not always possible.

`chem_scaffold_folds` first validates the original SMILES, removes `@`, `/`, and
`\` stereo annotations, then applies the native `murcko_scaffold` implementation
and groups equal scaffold strings. Fold identity is deliberately achiral, so
stereoisomers with otherwise identical supported structure stay together. All
acyclic molecules share one empty-scaffold group. Invalid molecules fail the
request. This inherits the
native parser, aromaticity, and canonicalization limits described in the README;
it does not promise grouping identical to another chemistry toolkit for every
molecule. Scaffold grouping reduces a particular source of structural leakage; it
does not remove duplicate-measurement, time, assay, target, or other leakage. Keep
related measurements together using an appropriate explicit grouping when needed.

## Experimental design and response surfaces

| SQL function / operation | Data rows | Parameters in order | Main JSON fields |
| --- | --- | --- | --- |
| `chem_factorial_design` / `factorial_design` | One `[low, high]` row per factor | `[center_points=0]` | `rows`, `coded_rows`, `run_count`, `factors` |
| `chem_central_composite_design` / `central_composite_design` | One `[low, high]` row per factor | `[alpha=2^(p/4), center_points=1]` | Same fields, plus `alpha` |
| `chem_d_optimal_design` / `d_optimal_design` | Candidate feature rows | `[n_runs, max_iterations=50, seed=0]` | `indices`, `rows`, `log_determinant`, `means`, `scales`, `exchanges`, `converged` |
| `chem_response_surface` / `response_surface` | `[features..., target]` | `[lambda=0]` | `terms`, `means`, `scales`, `coefficients`, `model_parameters`, `fitted`, `rmse`, `r2` |
| `chem_response_surface_predict` / `response_surface_predict` | Raw feature rows | Saved `model_parameters` | `predictions` |

Bounds require `low < high`. A full factorial emits all `2^p` low/high
combinations plus requested centers. The central composite design adds `2p`
axial runs and centers to that factorial. Coded coordinates map to original units
using `x=(low+high)/2 + z*(high-low)/2`. **For `alpha > 1`, axial runs extend beyond
the supplied factorial bounds.** These are not hard feasibility constraints.
Choose a feasible face-centered design (`alpha=1`) or supply explicitly feasible
candidates to D-optimal selection when bounds must be respected. Center counts
must be nonnegative integers, and `alpha` must be positive.

D-optimal selection adds an intercept to the supplied feature matrix, initializes
a full-rank subset, adds rows greedily, and exchanges selected and unselected rows
to increase `log(det(X'X))`. Exactly duplicate candidate feature rows are collapsed;
returned indices refer to the original input and contain no repeats. The selected
raw rows must come from the supplied candidate set, so constraints can be encoded
by filtering that set beforehand. Encode desired squares or interactions in the
candidate matrix yourself when optimizing information for that model.

The reported determinant uses mean-centered features divided by each feature's
half-range; `means` and `scales` document that basis. Translation and nonsingular
feature rescaling preserve which subset maximizes the determinant, while changing
its numeric value. Rank-deficient candidates or fewer than `p+1` requested runs
are rejected. This is a local exchange procedure, not a certificate of the global
D-optimum. `converged=false` means the exchange iteration budget was exhausted or
no iterations were requested. Output indices are sorted in original input order.

Response surfaces fit the full quadratic in normalized variables
`z[j]=(x[j]-means[j])/scales[j]`, using arithmetic means and half-ranges (constant
features use scale `1`). For `p` features there are
`q=1+2p+p(p-1)/2` terms, ordered:

1. Intercept.
2. `z0, z1, ...`.
3. `z0^2, z1^2, ...`.
4. `z0*z1, z0*z2, ...`, iterating pairs with `i<j` lexicographically.

Column-equilibrated, pivoted QR solves the model. A positive `lambda` adds a Ridge
penalty on all normalized non-intercept terms and can identify an otherwise
underdetermined fit. `lambda=0` requires numerical full rank. Coefficients remain
in this normalized quadratic basis, unlike the linear regression coefficients.
`model_parameters` contains `[p, means..., scales..., coefficients...]` and can be
passed unchanged to `chem_response_surface_predict` with raw feature rows.
`fitted`, RMSE, and R² describe the training observations; no validation split,
confidence intervals, stationary-point search, or constrained response optimizer
is implied. For a constant response, this operation reports R² `1` for an exact
fit and `0` otherwise.

## Candidate-pool optimization

| SQL function / operation | Data rows | Required parameters in order | Main JSON fields |
| --- | --- | --- | --- |
| `chem_gp_predict` / `gp_predict` | Training prefix `[features..., observed_target]`, then candidates `[features..., dummy_target]` | `[n_train, length_scale, signal_variance, noise_variance, minimize, xi]` | `predictions`, `variances`, `expected_improvement`, `ranked_candidate_indices`, `feature_means`, `feature_scales`, `target_mean` |
| `chem_acquisition` / `acquisition` | `[predictive_mean, predictive_sigma]` | `[best, minimize, xi]` | `expected_improvement` |
| `chem_diversity_select` / `diversity_select` | Existing feature prefix, then candidate feature rows | `[n_select, n_existing]` | `selected_indices` |

The Gaussian process uses an isotropic RBF kernel after feature centering and
population-standard-deviation scaling fitted only on the observed prefix:

```text
k(a,b) = signal_variance * exp(-0.5 * sum_j ((z[a,j]-z[b,j])/length_scale)^2)
K_observed = k(X_train,X_train) + noise_variance * I
```

`length_scale` and `signal_variance` must be positive; `noise_variance` is a
nonnegative **variance**, not a standard deviation. A Cholesky solve rejects
numerically singular covariance matrices; repeated training inputs may need
positive noise. Target centering uses the observed mean; targets are not scaled.
Thus signal and noise variances are in squared target units. Hyperparameters are
supplied explicitly; there is no hidden likelihood optimization. Candidate dummy
targets must be finite (use `0.0`) but are ignored by fitting, incumbent selection,
and acquisition.

Predictions, latent posterior variances, and EI arrays describe only the candidate
suffix, in its original order. The variance excludes observation noise for a
future measurement. `ranked_candidate_indices` orders candidates by decreasing
expected improvement, with lower candidate index first on a tie. These indices
start at zero within the candidate suffix, not within the full training-plus-
candidate matrix. GP inference does not automatically remove candidate duplicates
or points already observed; filter them before the call if required.

Set `minimize=1` for minimization or `0` for maximization. `xi >= 0` is the
improvement margin. With `delta=best-mean-xi` for minimization or
`delta=mean-best-xi` for maximization, EI is
`delta*Phi(delta/sigma) + sigma*phi(delta/sigma)`, and `max(delta,0)` if `sigma=0`.
`chem_acquisition` takes predictive **standard deviation** `sigma`, not variance.
Its normal CDF uses convergent incomplete-gamma series and continued fractions
to preserve small tail probabilities. Suggestions search the supplied finite
pool; they do not optimize a continuous domain or run experiments.

Diversity selection is deterministic farthest-first selection under squared
Euclidean distances in the **raw supplied feature units**, without GP scaling or
an acquisition score. It excludes the existing prefix and exact repeated
positions. If no existing rows are supplied, it starts with the first candidate.
It fails if too few distinct new positions remain. Returned indices are zero
based in the **full matrix**, unlike the GP ranking's candidate-relative indices.
Scale or weight heterogeneous descriptors before this distance calculation.

## Statistical methods

| SQL function / operation | Data rows | Required parameters | Results and scope |
| --- | --- | --- | --- |
| `chem_summarize` / `summarize` | Numeric feature rows | `[]` | Column `means`, `sample_variances`, `standard_deviations`, `minima`, `maxima`; sample variance/SD are `null` for one row |
| `chem_welch_t_test` / `welch_t_test` | `[value, group]`, group `0` or `1` | `[]` | Two-sided Welch test: group sizes/means, group-0 minus group-1 `mean_difference`, `t_statistic`, Satterthwaite `degrees_of_freedom`, `p_value` |
| `chem_one_way_anova` / `one_way_anova` | `[value, group_id]` | `[]` | One-way between-groups F test: sums of squares, degrees of freedom, `f_statistic`, `p_value` |
| `chem_adjust_pvalues` / `adjust_pvalues` | `[p_value]` in `[0,1]` | `[method]`: `0`=Holm, `1`=Benjamini-Hochberg | `adjusted_p_values` in original input order |
| `chem_process_capability` / `process_capability` | `[value]` | `[lower_specification_limit, upper_specification_limit, within_process_sigma (optional)]` | Mean, overall sample SD, `pp`/`ppk`; `cp`/`cpk` only with supplied within-process sigma; descriptive mean ± 3 SD reference limits |
| `chem_kaplan_meier` / `kaplan_meier` | `[time, event]`; time nonnegative, event `1`=event, `0`=right-censored | `[]` | `times`, `at_risk`, `events`, `censored`, `survival`, `greenwood_variance`, `median_survival` |
| `chem_random_intercept` / `random_intercept` | `[value, group_id]` | `[]` | Balanced one-way random-intercept ANOVA variance estimator: `within_group_variance`, `between_group_variance`, `icc`, grand mean and group counts |

Summaries use sample variance `s²=sum_i (x_i-mean)²/(n-1)`. For Welch, let
`v_g=s_g²/n_g`; then `t=(mean_0-mean_1)/sqrt(v_0+v_1)` and the degrees of freedom
are `(v_0+v_1)² / (v_0²/(n_0-1)+v_1²/(n_1-1))`. One-way ANOVA uses
`SS_between=sum_g n_g*(mean_g-grand_mean)²`, `SS_within=sum_g sum_i (x_gi-mean_g)²`,
and `F=(SS_between/(k-1))/(SS_within/(n-k))`.

Welch requires at least two observations per group and fails when both group
variances are zero. ANOVA needs at least two groups, positive residual degrees of
freedom, and nonzero within-group variation. Its conventional F test assumes an
appropriate independent-error, common-variance model; it is not Welch ANOVA.
Welch/ANOVA p-values use numerical regularized incomplete beta calculations. No
paired tests, nonparametric tests, factorial ANOVA, Tukey comparisons, or automatic
assumption selection are implemented.

Capability computes `pp=(USL-LSL)/(6*s)` and
`ppk=min(USL-mean,mean-LSL)/(3*s)` using the **overall sample SD** `s`. At least
two observations and a positive overall SD are required. Without the optional
third parameter, `cp`/`cpk` and `within_process_sigma` are JSON `null`. When a
positive within-process sigma is supplied, Cp/Cpk use that value instead of the
overall SD; the function does not estimate within-process sigma from subgroups.
`individual_lower_reference_limit` and `individual_upper_reference_limit` are
descriptive `mean ± 3*s`, not validated control-chart limits. This does not assess
process stability, normality, autocorrelation, subgroup structure, or
measurement-system adequacy.

Kaplan-Meier handles right censoring, with `S(t)=product_(u<=t)(1-d_u/r_u)`
for event counts `d_u` and at-risk counts `r_u`. Greenwood variance is
`S(t)² * sum_(u<=t) d_u/(r_u*(r_u-d_u))`, with variance `0` after survival reaches
zero. At tied times, the risk set includes both events and same-time censors
before the event decrement. The median is the first
time survival reaches `0.5`, or JSON `null` if it never does. Greenwood variance
is returned without confidence intervals. There is no left truncation, interval
censoring, covariate model, log-rank test, Cox model, or general lifetime-model
fitter.

Random-intercept estimation requires at least two balanced groups with at least
two replicates each. Residual variance is `MS_within`; the group variance is
`max((MS_between-MS_within)/replicates,0)`. ICC is group variance divided by total
variance, or `null` when total variance is zero. This is not a general mixed-model
solver: no unbalanced groups, fixed covariates, random slopes, REML fitting,
BLUPs, or variance-component confidence intervals are supplied.

## Resource limits

The SQL/C ABI boundary accepts at most 100,000 rows, 1,024 columns, 1,000,000
numeric cells, and 16,384 parameter values, and limits output to 64 MiB. The
following stricter operation limits also apply; finite inputs can still fail if
intermediate calculations overflow or numerical rank is inadequate.

| Operations | Additional limits |
| --- | --- |
| PCA / k-means | At most 256 features; absolute input magnitude <= `1e150`. PCA `n*p*p <= 500,000,000`; k-means `n*p*k*max_iterations <= 500,000,000`. |
| Ridge / PLS1 / CV | At most 256 predictor features; fit `n*p*p <= 100,000,000`; CV `n_folds*n*p*p <= 100,000,000`. |
| Scaffold folds | At most 100,000 molecules, each at most 4,096 SMILES bytes, and at most 4 MiB of SMILES per call; native scaffold serialization is bounded to 65,536 bytes. |
| DOE / response surface | At most 16 factors/features; matrix at most 65,536 rows and 32 columns. Designs at most 65,536 runs and 1,000,000 cells. |
| D-optimal | At most 2,048 candidates, 256 runs, 1,000 exchange iterations; `candidates*runs*(p+1)^2*max(iterations,1) <= 500,000,000`. |
| Response surface | With `q` quadratic terms, `n*q*q <= 100,000,000` and `(n+q)*q <= 1,000,000`. |
| GP / acquisition / diversity | Matrix at most 130 columns. GP: `1..512` observed rows followed by `1..10000` candidates. Diversity: at most 10,000 selections and `(n_select+n_existing)*n*p <= 20,000,000`. |
| Statistical methods | General SQL row/column/cell bounds; additional group/variance requirements above. |

Seeds, group IDs, and fold IDs must be exact nonnegative integers no larger than
`9007199254740991`. Numeric matrices are materialized for one call; filter large
source tables and select useful features before building nested lists. These APIs
are not streaming or distributed solvers.

## End-to-end SQL example

The following script uses invented target values only to show the API flow. Run
it in one connection after loading a build containing these functions. The JSON
extraction examples also require DuckDB's `json` extension (`INSTALL json; LOAD
json;` if it is not built in or already available). All aggregations have explicit
row order. The example uses molecule descriptors for
exploration/prediction, scaffold-grouped validation, a separate feasible process
DOE, and candidate-pool GP suggestions.

```sql
-- LOAD ducksmiles; -- or load your locally built extension first
-- INSTALL json; LOAD json; -- if JSON functions are not already available
CREATE TEMP TABLE ds_observed AS
SELECT * FROM (VALUES
  (1, 'c1ccccc1', 1.2),
  (2, 'Cc1ccccc1', 1.4),
  (3, 'C1CCCCC1', 0.7),
  (4, 'CC1CCCCC1', 0.8),
  (5, 'c1ccncc1', 2.0),
  (6, 'Cc1ccncc1', 2.2),
  (7, 'CCO', 0.3),
  (8, 'CCCO', 0.5)
) AS t(id, smiles, target);

CREATE TEMP TABLE ds_features AS
SELECT id, smiles, target::DOUBLE AS target,
       mol_weight(smiles) AS mw,
       logp_crippen(smiles) AS logp,
       tpsa(smiles) AS polar
FROM ds_observed;

CREATE TEMP TABLE ds_pca AS
SELECT chem_pca(list([mw,logp,polar] ORDER BY id), [2.0,1.0]) AS model
FROM ds_features;
SELECT json_extract(model, '$.explained_variance_ratio') AS variance_ratio
FROM ds_pca;

CREATE TEMP TABLE ds_ridge AS
SELECT chem_ridge(list([mw,logp,polar,target] ORDER BY id), [1.0]) AS model
FROM ds_features;
SELECT json_extract(model, '$.metrics') AS training_metrics FROM ds_ridge;

CREATE TEMP TABLE ds_folds AS
SELECT chem_scaffold_folds(list(smiles ORDER BY id), [3.0,42.0]) AS model
FROM ds_features;
CREATE TEMP TABLE ds_validation AS
SELECT f.*,
       (json_extract(g.model, '$.folds')::DOUBLE[])
         [row_number() OVER (ORDER BY id)] AS fold
FROM ds_features f CROSS JOIN ds_folds g;
SELECT chem_cross_validate(
         list([mw,logp,polar,target,fold] ORDER BY id), [0.0,1.0]) AS ridge_cv,
       chem_cross_validate(
         list([mw,logp,polar,target,fold] ORDER BY id), [1.0,2.0]) AS pls_cv
FROM ds_validation;

-- Temperature and composition bounds; alpha=1 keeps the axial points feasible.
SELECT chem_factorial_design([[20.0,80.0],[0.1,0.9]], [1.0]) AS factorial;
CREATE TEMP TABLE ds_ccd AS
SELECT chem_central_composite_design(
         [[20.0,80.0],[0.1,0.9]], [1.0,1.0]) AS design;
SELECT json_extract(design, '$.rows') AS process_runs FROM ds_ccd;
SELECT chem_d_optimal_design(
         [[20.0,0.1],[20.0,0.9],[80.0,0.1],[80.0,0.9],[50.0,0.5]],
         [4.0,30.0,42.0]) AS selected_feasible_runs;

-- Synthetic process measurements demonstrate a reusable full quadratic model.
CREATE TEMP TABLE ds_process AS
SELECT run_id, point,
       3.0 + 0.1*point[1] - 2.0*point[2]
           + 0.01*point[1]*point[1] + point[2]*point[2]
           + 0.02*point[1]*point[2] AS response
FROM (
  SELECT i AS run_id, points[i] AS point
  FROM (SELECT json_extract(design, '$.rows')::DOUBLE[][] AS points FROM ds_ccd)
  CROSS JOIN range(1, len(points)+1) AS t(i)
);
CREATE TEMP TABLE ds_surface AS
SELECT chem_response_surface(
         list(list_append(point,response) ORDER BY run_id), [0.0]) AS model
FROM ds_process;
SELECT chem_response_surface_predict(
         [[40.0,0.3],[60.0,0.7]],
         json_extract(model, '$.model_parameters')::DOUBLE[]) AS predictions
FROM ds_surface;

CREATE TEMP TABLE ds_candidates AS
SELECT id, smiles, mol_weight(smiles) AS mw,
       logp_crippen(smiles) AS logp, tpsa(smiles) AS polar
FROM (VALUES (101,'CCc1ccccc1'), (102,'CC1CCCC1'),
             (103,'Oc1ccncc1'), (104,'CCCCO')) AS t(id,smiles);
WITH x AS (
  SELECT list([mw,logp,polar] ORDER BY id) AS features FROM ds_candidates
)
SELECT chem_predict(x.features,
         json_extract(m.model, '$.stable_predict_parameters')::DOUBLE[])
       AS ridge_candidate_predictions
FROM x CROSS JOIN ds_ridge m;

-- Prefix training rows, then candidate rows with ignored finite dummy targets.
CREATE TEMP TABLE ds_gp_rows AS
SELECT 0 AS phase, id, mw, logp, polar, target FROM ds_features
UNION ALL
SELECT 1 AS phase, id, mw, logp, polar, 0.0 AS target FROM ds_candidates;
CREATE TEMP TABLE ds_gp AS
SELECT chem_gp_predict(
         list([mw,logp,polar,target] ORDER BY phase,id),
         [8.0,1.0,1.0,0.01,0.0,0.0]) AS model
FROM ds_gp_rows;
WITH ranking AS (
  SELECT unnest(json_extract(model, '$.ranked_candidate_indices')::INTEGER[])
         AS candidate_index FROM ds_gp
), indexed_candidates AS (
  SELECT *, row_number() OVER (ORDER BY id)-1 AS candidate_index
  FROM ds_candidates
)
SELECT c.id, c.smiles,
       (json_extract(g.model, '$.predictions')::DOUBLE[])
         [r.candidate_index+1] AS predicted_target,
       (json_extract(g.model, '$.expected_improvement')::DOUBLE[])
         [r.candidate_index+1] AS expected_improvement
FROM ranking r JOIN indexed_candidates c USING (candidate_index)
CROSS JOIN ds_gp g
ORDER BY expected_improvement DESC, c.id;
```

## Verification and extension boundaries

The Rust crate is `crates/analytics`; its tests cover independently calculated
numerical cases, rank/shape validation, fold-local preprocessing, design
properties, candidate indexing, and C buffer ownership. Use one build directory:

```sh
CARGO_TARGET_DIR=/private/tmp/rusta-ducksmiles-analytics-target cargo test --workspace
CARGO_TARGET_DIR=/private/tmp/rusta-ducksmiles-analytics-target cargo clippy -p ducksmiles_analytics --all-targets -- -D warnings
```

The SQL wrappers register typed functions and preserve nested-list validity and
DuckDB selection vectors. The C ABI receives a row-major matrix, returns an owned
UTF-8 JSON buffer, and requires freeing it with `ds_analytics_free` using its exact
length. Rust calls expose error text as `Result`; the SQL wrapper deliberately
returns `NULL` on failure. A successful Rust test does not replace loading the
matching compiled DuckDB extension and exercising its SQL boundary.
