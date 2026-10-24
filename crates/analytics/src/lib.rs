//! Dependency-free numerical analytics for the duckSMILES SQL extension.
pub mod advanced_design;
pub mod arima;
pub mod causal;
pub mod censored_regression;
pub mod degradation;
pub mod design;
pub mod design_classical;
pub mod distributions;
pub mod dynamic_series;
pub mod evaluation;
pub mod explore;
pub mod functional_preprocess;
pub mod genomics;
pub mod genreg_ordinal;
pub mod glm;
pub mod glmm;
pub mod inference;
pub mod irt;
pub mod kinship;
pub mod latent;
pub mod latent_inference;
pub mod learning;
pub mod market;
pub mod mixed;
pub mod mixed_advanced;
pub mod multivariate;
pub mod nonlinear;
pub mod omics;
pub mod optimize;
pub mod pls_multivariate;
pub mod quality;
pub mod regression;
pub mod reliability;
pub mod resampling;
pub mod robust_learning;
pub mod sem;
pub mod spatial;
pub mod spatial_estimation;
pub mod spectral;
pub mod state_space;
pub mod static_reports;
pub mod statistics;
pub mod survival;
pub mod text_analysis;
pub mod timeseries;
pub mod zero_inflation;

/// # Safety
/// Pointers must reference the indicated readable lengths; `output_len` is writable.
#[no_mangle]
pub unsafe extern "C" fn ds_analytics_run(
    operation: *const u8,
    operation_len: usize,
    data: *const f64,
    rows: usize,
    columns: usize,
    parameters: *const f64,
    parameters_len: usize,
    output_len: *mut usize,
) -> *mut u8 {
    use std::{
        panic::{catch_unwind, AssertUnwindSafe},
        ptr, slice,
    };
    if output_len.is_null() {
        return ptr::null_mut();
    }
    *output_len = 0;
    let count = match rows.checked_mul(columns) {
        Some(count) if count > 0 && count <= 1_000_000 && rows <= 100_000 && columns <= 1024 => {
            count
        }
        _ => return ptr::null_mut(),
    };
    if operation.is_null()
        || operation_len == 0
        || operation_len > 64
        || data.is_null()
        || parameters_len > 16384
        || (parameters_len > 0 && parameters.is_null())
    {
        return ptr::null_mut();
    }
    let outcome = catch_unwind(AssertUnwindSafe(|| -> Result<Vec<u8>, String> {
        // Lengths are checked before constructing slices. Caller guarantees allocation validity.
        let operation = std::str::from_utf8(slice::from_raw_parts(operation, operation_len))
            .map_err(|_| "operation is not UTF-8".to_string())?;
        let flat = slice::from_raw_parts(data, count);
        let data: Vec<Vec<f64>> = flat.chunks_exact(columns).map(|row| row.to_vec()).collect();
        let parameters = if parameters_len == 0 {
            &[]
        } else {
            slice::from_raw_parts(parameters, parameters_len)
        };
        let result = execute(operation, &data, parameters)?;
        if result.len() > 64 * 1024 * 1024 {
            return Err("analytics result exceeds 64 MiB".into());
        }
        Ok(result.into_bytes())
    }));
    match outcome {
        Ok(Ok(result)) => {
            let result = result.into_boxed_slice();
            *output_len = result.len();
            Box::into_raw(result) as *mut u8
        }
        _ => ptr::null_mut(),
    }
}

/// # Safety
/// The pointer must be a live buffer returned by `ds_analytics_run`, with its exact length.
#[no_mangle]
pub unsafe extern "C" fn ds_analytics_free(data: *mut u8, len: usize) {
    if !data.is_null() {
        drop(Box::from_raw(std::ptr::slice_from_raw_parts_mut(data, len)));
    }
}

pub fn execute(operation: &str, data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    match operation {
        "pls2_fit" | "pls2_predict" | "plsda_fit" | "plsda_predict" | "pls2_validate"
        | "plsda_validate" => pls_multivariate::run(operation, data, parameters),
        "zi_binomial"
        | "zi_beta_binomial"
        | "zi_poisson"
        | "zi_negative_binomial"
        | "zi_gamma"
        | "zi_predict" => zero_inflation::run(operation, data, parameters),
        "genreg_ordinal_fit" | "genreg_ordinal_predict" => {
            genreg_ordinal::run(operation, data, parameters)
        }
        "report_histogram_svg" | "report_box_plot_svg" => {
            static_reports::run(operation, data, parameters)
        }
        "censored_normal"
        | "censored_exponential"
        | "censored_gamma"
        | "censored_weibull"
        | "censored_lognormal"
        | "censored_beta"
        | "censored_predict" => censored_regression::run(operation, data, parameters),
        "spatial_fit" | "spatial_fit_predict" => {
            spatial_estimation::run(operation, data, parameters)
        }
        "kalman_smoother" | "state_space_smoothing" | "state_space_predict" => {
            state_space::run(operation, data, parameters)
        }
        "arima_forecast_intervals"
        | "cross_correlation"
        | "arima_prewhitening"
        | "transfer_function_model"
        | "remove_linear_trend"
        | "remove_cycle" => dynamic_series::run(operation, data, parameters),
        "functional_grid"
        | "functional_bin"
        | "functional_thin"
        | "functional_snip"
        | "functional_arwls"
        | "functional_baseline"
        | "functional_dtw_matrix"
        | "functional_golden_match" => functional_preprocess::run(operation, data, parameters),
        "text_dtm"
        | "text_weight"
        | "text_weight_predict"
        | "text_lsa"
        | "text_projection"
        | "text_topic_svd"
        | "text_latent_class"
        | "text_latent_class_predict" => text_analysis::run(operation, data, parameters),
        "spatial_covariance" | "spatial_variogram" | "spatial_gls" | "spatial_predict" => {
            spatial::run(operation, data, parameters)
        }
        "kinship_model" | "kinship_predict" => kinship::run(operation, data, parameters),
        "irt_fit" | "irt_predict" | "irt_score" => irt::run(operation, data, parameters),
        "cfa_model" | "cfa_scores" | "sem_fit_indices" | "sem_ram" => {
            sem::run(operation, data, parameters)
        }
        "degradation_path"
        | "degradation_predict"
        | "stability_analysis"
        | "destructive_degradation"
        | "destructive_predict"
        | "cumulative_damage_fit"
        | "cumulative_damage_predict" => degradation::run(operation, data, parameters),
        "rcbd_design"
        | "balanced_incomplete_block"
        | "resolvable_incomplete_block"
        | "split_plot_design"
        | "strip_plot_design"
        | "taguchi_orthogonal_array" => design_classical::run(operation, data, parameters),
        "association_rules" | "basket_affinity" | "item_reliability" => {
            market::run(operation, data, parameters)
        }
        "propensity_fit" | "propensity_predict" | "causal_balance" | "propensity_match"
        | "ipw_effect" | "aipw_effect" | "uplift_curve" => causal::run(operation, data, parameters),
        "assay_normalize"
        | "quantile_normalize"
        | "quantile_normalize_predict"
        | "marker_admixture"
        | "admixture_predict" => omics::run(operation, data, parameters),
        "glmm_binomial" | "glmm_poisson" | "glmm_predict" => glmm::run(operation, data, parameters),
        "repeated_covariance"
        | "repeated_model"
        | "random_slope_model"
        | "random_slope_predict" => mixed_advanced::run(operation, data, parameters),
        "aft_exponential" | "aft_weibull" | "aft_lognormal" | "aft_predict" | "nhpp_powerlaw"
        | "nhpp_predict" | "reliability_blocks" => reliability::run(operation, data, parameters),
        "marker_statistics"
        | "marker_impute"
        | "marker_impute_predict"
        | "marker_relatedness"
        | "distance_matrix"
        | "marker_simulate"
        | "marker_cross"
        | "genomic_blup"
        | "genomic_blup_predict" => genomics::run(operation, data, parameters),
        "svm_fit"
        | "svm_predict"
        | "quantile_regression"
        | "huber_regression"
        | "robust_predict" => robust_learning::run(operation, data, parameters),
        "manova"
        | "canonical_correlation_inference"
        | "multiple_correspondence_analysis"
        | "factor_analysis_ml" => latent_inference::run(operation, data, parameters),
        "spectral_snv"
        | "spectral_msc"
        | "savitzky_golay"
        | "curve_derivative"
        | "spectral_detrend"
        | "haar_transform"
        | "fourier_transform"
        | "dynamic_time_warping"
        | "functional_pca"
        | "peak_find" => spectral::run(operation, data, parameters),
        "curve_fit" | "curve_predict" | "bspline_basis" | "spline_fit" | "spline_predict"
        | "loess" | "loess_predict" => nonlinear::run(operation, data, parameters),
        "regression_metrics"
        | "confusion_matrix"
        | "binary_metrics"
        | "roc_curve"
        | "precision_recall_curve"
        | "lift_curve"
        | "calibration_curve"
        | "profit_decision"
        | "ensemble_predict"
        | "model_compare"
        | "holdout_split"
        | "rolling_origin" => evaluation::run(operation, data, parameters),
        "mds"
        | "correspondence_analysis"
        | "factor_analysis"
        | "canonical_correlation"
        | "gaussian_mixture"
        | "gaussian_mixture_predict"
        | "self_organizing_map"
        | "som_predict" => latent::run(operation, data, parameters),
        "box_behnken"
        | "latin_hypercube"
        | "mixture_lattice"
        | "fractional_factorial"
        | "orthogonal_screening"
        | "randomize_design"
        | "design_evaluate" => advanced_design::run(operation, data, parameters),
        "arima" => arima::run(operation, data, parameters),
        "mixed_model" | "mixed_predict" => mixed::run(operation, data, parameters),
        "cox_ph" | "cox_predict" | "relative_risk_predict" => {
            survival::run(operation, data, parameters)
        }
        "decision_tree" | "random_forest" | "gradient_boosting" | "mlp" => {
            learning::run(operation, data, parameters)
        }
        "acf"
        | "pacf"
        | "difference"
        | "moving_average"
        | "exponential_smoothing"
        | "holt"
        | "ar_model"
        | "periodogram"
        | "seasonal_decompose"
        | "ljung_box" => timeseries::run(operation, data, parameters),
        "covariance_matrix"
        | "correlation_matrix"
        | "hierarchical_cluster"
        | "mahalanobis"
        | "discriminant_analysis"
        | "gaussian_naive_bayes"
        | "knn_predict"
        | "dbscan" => multivariate::run(operation, data, parameters),
        "fit_normal"
        | "fit_lognormal"
        | "fit_exponential"
        | "fit_poisson"
        | "fit_gamma"
        | "fit_weibull"
        | "normality_test"
        | "kde_gaussian"
        | "box_cox"
        | "survival_exponential"
        | "survival_weibull"
        | "logrank_test" => distributions::run(operation, data, parameters),
        "bootstrap_mean"
        | "permutation_test"
        | "sample_rows"
        | "simulate_normal"
        | "power_z"
        | "sample_size_mean"
        | "sample_size_proportion"
        | "reliability_demonstration"
        | "odds_ratio"
        | "relative_risk"
        | "cohen_d" => resampling::run(operation, data, parameters),
        "one_sample_t"
        | "paired_t"
        | "pooled_t"
        | "variance_test"
        | "variance_ratio_test"
        | "proportion_test"
        | "contingency_test"
        | "fisher_exact"
        | "mann_whitney"
        | "wilcoxon_signed_rank"
        | "kruskal_wallis"
        | "levene_test"
        | "equivalence_t" => inference::run(operation, data, parameters),
        "linear_model"
        | "logistic"
        | "poisson"
        | "glm_predict"
        | "elastic_net"
        | "polynomial_features" => glm::run(operation, data, parameters),
        "ir_chart"
        | "xbar_r_chart"
        | "xbar_s_chart"
        | "p_chart"
        | "np_chart"
        | "c_chart"
        | "u_chart"
        | "ewma_chart"
        | "cusum_chart"
        | "moving_average_chart"
        | "western_electric"
        | "gage_rr"
        | "type1_gage"
        | "agreement_kappa"
        | "acceptance_oc"
        | "levey_jennings"
        | "laney_p_chart"
        | "laney_u_chart" => quality::run(operation, data, parameters),
        "pca" | "kmeans" => explore::run(operation, data, parameters),
        "ridge" | "pls" | "predict" | "cross_validate" | "group_folds" => {
            regression::run(operation, data, parameters)
        }
        "factorial_design"
        | "central_composite_design"
        | "d_optimal_design"
        | "response_surface"
        | "response_surface_predict" => design::run(operation, data, parameters),
        "summarize" | "welch_t_test" | "one_way_anova" | "adjust_pvalues"
        | "process_capability" | "kaplan_meier" | "random_intercept" => {
            statistics::run(operation, data, parameters)
        }
        _ => optimize::run(operation, data, parameters),
    }
}
