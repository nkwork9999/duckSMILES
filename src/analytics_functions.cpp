#include "analytics.h"
#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#include "ducksmiles.h"
#include "ducksmiles_compat.hpp"
#include "ducksmiles_extension.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

namespace duckdb {
namespace {
struct AnalyticsBuffer {
  uint8_t *data;
  size_t size;
  ~AnalyticsBuffer() { ds_analytics_free(data, size); }
};

bool ReadParameters(const Value &value, std::vector<double> &parameters) {
  if (value.IsNull())
    return false;
  const auto &children = ListValue::GetChildren(value);
  if (children.size() > 16384)
    return false;
  for (const auto &child : children) {
    if (child.IsNull())
      return false;
    double number = child.GetValue<double>();
    if (!std::isfinite(number))
      return false;
    parameters.push_back(number);
  }
  return true;
}

bool ReadMatrix(const Value &value, std::vector<double> &data, size_t &rows,
                size_t &columns) {
  if (value.IsNull())
    return false;
  const auto &outer = ListValue::GetChildren(value);
  rows = outer.size();
  if (rows == 0 || rows > 100000 || outer[0].IsNull())
    return false;
  columns = ListValue::GetChildren(outer[0]).size();
  if (columns == 0 || columns > 1024 || rows > 1000000 / columns)
    return false;
  data.reserve(rows * columns);
  for (const auto &row : outer) {
    if (row.IsNull())
      return false;
    const auto &children = ListValue::GetChildren(row);
    if (children.size() != columns)
      return false;
    for (const auto &child : children) {
      if (child.IsNull())
        return false;
      double number = child.GetValue<double>();
      if (!std::isfinite(number))
        return false;
      data.push_back(number);
    }
  }
  return true;
}

// All rowwise access goes through DuckDB Value, preserving dictionary/constant
// selection and nested-list validity. Rust owns the result until copied here.
void RunAnalytics(DataChunk &args, Vector &result,
                  const char *fixed_operation) {
  result.SetVectorType(VectorType::FLAT_VECTOR);
  auto output = ducksmiles_compat::ResultData<string_t>(result, args.size());
  auto &validity = ducksmiles_compat::ResultValidity(result);
  for (idx_t row = 0; row < args.size(); row++) {
    const idx_t offset = fixed_operation ? 0 : 1;
    std::string operation;
    if (fixed_operation)
      operation = fixed_operation;
    else {
      auto value = args.GetValue(0, row);
      if (value.IsNull()) {
        validity.SetInvalid(row);
        continue;
      }
      operation = value.GetValue<std::string>();
    }
    std::vector<double> data, parameters;
    size_t rows = 0, columns = 0;
    if (!ReadMatrix(args.GetValue(offset, row), data, rows, columns) ||
        !ReadParameters(args.GetValue(offset + 1, row), parameters)) {
      validity.SetInvalid(row);
      continue;
    }
    size_t length = 0;
    auto bytes =
        ds_analytics_run(reinterpret_cast<const uint8_t *>(operation.data()),
                         operation.size(), data.data(), rows, columns,
                         parameters.data(), parameters.size(), &length);
    AnalyticsBuffer buffer{bytes, length};
    if (!bytes) {
      validity.SetInvalid(row);
      continue;
    }
    output[row] = StringVector::AddString(
        result, reinterpret_cast<const char *>(bytes), length);
  }
}

void Analyze(DataChunk &args, ExpressionState &, Vector &result) {
  RunAnalytics(args, result, nullptr);
}

#define ANALYTICS_CALLBACK(Name, Operation)                                    \
  void Name(DataChunk &args, ExpressionState &, Vector &result) {              \
    RunAnalytics(args, result, Operation);                                     \
  }
ANALYTICS_CALLBACK(Pca, "pca")
ANALYTICS_CALLBACK(Kmeans, "kmeans")
ANALYTICS_CALLBACK(Ridge, "ridge")
ANALYTICS_CALLBACK(Pls, "pls")
ANALYTICS_CALLBACK(Predict, "predict")
ANALYTICS_CALLBACK(CrossValidate, "cross_validate")
ANALYTICS_CALLBACK(GroupFolds, "group_folds")
ANALYTICS_CALLBACK(FactorialDesign, "factorial_design")
ANALYTICS_CALLBACK(CentralCompositeDesign, "central_composite_design")
ANALYTICS_CALLBACK(DOptimalDesign, "d_optimal_design")
ANALYTICS_CALLBACK(ResponseSurface, "response_surface")
ANALYTICS_CALLBACK(ResponseSurfacePredict, "response_surface_predict")
ANALYTICS_CALLBACK(GpPredict, "gp_predict")
ANALYTICS_CALLBACK(Acquisition, "acquisition")
ANALYTICS_CALLBACK(DiversitySelect, "diversity_select")
ANALYTICS_CALLBACK(Summarize, "summarize")
ANALYTICS_CALLBACK(WelchTest, "welch_t_test")
ANALYTICS_CALLBACK(Anova, "one_way_anova")
ANALYTICS_CALLBACK(AdjustPvalues, "adjust_pvalues")
ANALYTICS_CALLBACK(ProcessCapability, "process_capability")
ANALYTICS_CALLBACK(KaplanMeier, "kaplan_meier")
ANALYTICS_CALLBACK(RandomIntercept, "random_intercept")
#define ANALYTICS_METHOD(Symbol, Name, Operation, Description, Example)        \
  ANALYTICS_CALLBACK(Symbol, Operation)
#include "analytics_registry.inc"
#undef ANALYTICS_METHOD
#undef ANALYTICS_CALLBACK

void ScaffoldFolds(DataChunk &args, ExpressionState &, Vector &result) {
  result.SetVectorType(VectorType::FLAT_VECTOR);
  auto output = ducksmiles_compat::ResultData<string_t>(result, args.size());
  auto &validity = ducksmiles_compat::ResultValidity(result);
  for (idx_t row = 0; row < args.size(); row++) {
    auto list = args.GetValue(0, row);
    std::vector<double> parameters;
    if (list.IsNull() || !ReadParameters(args.GetValue(1, row), parameters)) {
      validity.SetInvalid(row);
      continue;
    }
    const auto &molecules = ListValue::GetChildren(list);
    if (molecules.empty() || molecules.size() > 100000) {
      validity.SetInvalid(row);
      continue;
    }
    std::vector<std::string> keys;
    std::map<std::string, size_t> groups;
    bool valid = true;
    size_t total_bytes = 0;
    for (const auto &molecule : molecules) {
      if (molecule.IsNull()) {
        valid = false;
        break;
      }
      auto smiles = molecule.GetValue<std::string>();
      total_bytes += smiles.size();
      // Bound chemistry work before calling its existing parser/serializer.
      if (smiles.size() > 4096 || total_bytes > 4 * 1024 * 1024) {
        valid = false;
        break;
      }
      auto ptr = reinterpret_cast<const uint8_t *>(smiles.data());
      if (ds_mol_is_valid(ptr, smiles.size()) != 1) {
        valid = false;
        break;
      }
      // Fold identity is achiral. Remove raw stereo annotations before native
      // canonicalization, whose traversal order can itself depend on stereo.
      smiles.erase(std::remove_if(smiles.begin(), smiles.end(),
                                  [](char c) {
                                    return c == '@' || c == '/' || c == '\\';
                                  }),
                   smiles.end());
      ptr = reinterpret_cast<const uint8_t *>(smiles.data());
      int32_t needed = ds_murcko_scaffold(ptr, smiles.size(), nullptr, 0);
      if (needed < 0 || needed > 65536) {
        valid = false;
        break;
      }
      std::vector<uint8_t> buffer(static_cast<size_t>(needed) + 1);
      int32_t written =
          ds_murcko_scaffold(ptr, smiles.size(), buffer.data(), buffer.size());
      if (written != needed) {
        valid = false;
        break;
      }
      keys.emplace_back(reinterpret_cast<const char *>(buffer.data()),
                        static_cast<size_t>(written));
      groups.emplace(keys.back(), 0);
    }
    if (!valid) {
      validity.SetInvalid(row);
      continue;
    }
    size_t index = 0;
    for (auto &group : groups)
      group.second = index++;
    std::vector<double> data;
    for (const auto &key : keys)
      data.push_back(static_cast<double>(groups[key]));
    const std::string operation = "group_folds";
    size_t length = 0;
    auto bytes =
        ds_analytics_run(reinterpret_cast<const uint8_t *>(operation.data()),
                         operation.size(), data.data(), data.size(), 1,
                         parameters.data(), parameters.size(), &length);
    AnalyticsBuffer buffer{bytes, length};
    if (!bytes) {
      validity.SetInvalid(row);
      continue;
    }
    output[row] = StringVector::AddString(
        result, reinterpret_cast<const char *>(bytes), length);
  }
}

void Register(ExtensionLoader &loader, const char *name,
              scalar_function_t callback, const string &description,
              const string &example, bool generic = false,
              bool scaffold = false) {
  auto matrix = LogicalType::LIST(LogicalType::LIST(LogicalType::DOUBLE));
  auto numbers = LogicalType::LIST(LogicalType::DOUBLE);
  vector<LogicalType> arguments =
      generic ? vector<LogicalType>{LogicalType::VARCHAR, matrix, numbers}
              : vector<LogicalType>{
                    scaffold ? LogicalType::LIST(LogicalType::VARCHAR) : matrix,
                    numbers};
  ScalarFunction function(name, arguments, LogicalType::VARCHAR, callback);
  function.SetFallible();
  FunctionDescription documentation;
  documentation.parameter_types = arguments;
  documentation.parameter_names =
      generic ? vector<string>{"operation", "data", "parameters"}
              : vector<string>{scaffold ? "smiles" : "data", "parameters"};
  ducksmiles_compat::SetArgumentNames(function, documentation.parameter_names);
  documentation.description =
      description + " Returns versioned JSON; NULL for invalid, nonfinite, "
                    "rank-deficient or resource-limited inputs. See "
                    "ANALYTICS.md for parameter order and method limits.";
  documentation.examples = {example};
  documentation.categories = {"ducksmiles", "analytics"};
  CreateScalarFunctionInfo info(std::move(function));
  info.on_conflict = OnCreateConflict::ALTER_ON_CONFLICT;
  info.descriptions.push_back(std::move(documentation));
  loader.RegisterFunction(std::move(info));
}
} // namespace

void RegisterAnalyticsFunctions(ExtensionLoader &loader) {
#define ANALYTICS_METHOD(Symbol, Name, Operation, Description, Example)        \
  Register(loader, Name, Symbol, Description, Example);
#include "analytics_registry.inc"
#undef ANALYTICS_METHOD
  Register(loader, "chem_analyze", Analyze,
           "Run a named numerical analytics operation.",
           "chem_analyze('acquisition', [[2.0,1.0]], [2.0,0.0,0.0])", true);
  Register(loader, "chem_pca", Pca,
           "PCA: parameters [components,standardize 0/1].",
           "chem_pca([[1.0,2.0],[2.0,4.0],[3.0,6.0]], [1.0,0.0])");
  Register(loader, "chem_kmeans", Kmeans,
           "Deterministic k-means: [clusters,max_iterations,seed].",
           "chem_kmeans([[0.0],[1.0],[9.0],[10.0]], [2.0,100.0,42.0])");
  Register(loader, "chem_ridge", Ridge,
           "Ridge regression: rows [features...,target], parameters [lambda].",
           "chem_ridge([[0.0,1.0],[1.0,3.0],[2.0,5.0]], [1.0])");
  Register(
      loader, "chem_pls", Pls,
      "PLS1 regression: rows [features...,target], parameters [components].",
      "chem_pls([[0.0,1.0],[1.0,3.0],[2.0,5.0]], [1.0])");
  Register(loader, "chem_predict", Predict,
           "Linear prediction: saved stable_predict_parameters or "
           "[intercept,coefficients...].",
           "chem_predict([[1.0],[2.0]], [1.0,2.0])");
  Register(loader, "chem_cross_validate", CrossValidate,
           "Fold-local CV: rows [features...,target,fold], [method 0 ridge/1 "
           "PLS,hyperparameter].",
           "chem_cross_validate([[0.0,1.0,0.0],[1.0,3.0,0.0],[2.0,5.0,1.0],[3."
           "0,7.0,1.0]], [0.0,1.0])");
  Register(loader, "chem_group_folds", GroupFolds,
           "Balanced indivisible group folds: rows [group_id], parameters "
           "[n_folds,seed].",
           "chem_group_folds([[1.0],[1.0],[2.0],[3.0]], [2.0,42.0])");
  Register(loader, "chem_scaffold_folds", ScaffoldFolds,
           "Achiral Murcko scaffold grouping and balanced fold assignment. All "
           "acyclic molecules share one empty-scaffold group; grouping "
           "inherits native chemistry canonicalization limits.",
           "chem_scaffold_folds(['c1ccccc1','Cc1ccccc1','C1CCCCC1','CCO'], "
           "[2.0,42.0])",
           false, true);
  Register(
      loader, "chem_factorial_design", FactorialDesign,
      "Full two-level factorial: rows [low,high], parameters [center_points].",
      "chem_factorial_design([[10.0,20.0],[1.0,3.0]], [1.0])");
  Register(loader, "chem_central_composite_design", CentralCompositeDesign,
           "Central composite design: rows [low,high], parameters "
           "[alpha,center_points]. Axial points may exceed bounds.",
           "chem_central_composite_design([[10.0,20.0],[1.0,3.0]], [1.0,1.0])");
  Register(
      loader, "chem_d_optimal_design", DOptimalDesign,
      "Local D-optimal exchange over supplied candidate rows, with intercept: "
      "[n_runs,max_iterations,seed].",
      "chem_d_optimal_design([[-1.0],[-0.5],[0.5],[1.0]], [2.0,20.0,42.0])");
  Register(loader, "chem_response_surface", ResponseSurface,
           "Full quadratic response fit: rows [features...,target], parameters "
           "[lambda].",
           "chem_response_surface([[-1.0,2.0],[0.0,1.0],[1.0,4.0]], [0.0])");
  Register(loader, "chem_response_surface_predict", ResponseSurfacePredict,
           "Quadratic prediction using flattened model_parameters from "
           "chem_response_surface.",
           "chem_response_surface_predict([[0.0],[1.0]], "
           "[1.0,0.0,1.0,1.0,1.0,2.0])");
  Register(
      loader, "chem_gp_predict", GpPredict,
      "RBF GP and EI on candidate pool: training prefix then candidates; "
      "[n_train,length_scale,signal_variance,noise_variance,minimize,xi].",
      "chem_gp_predict([[0.0,2.0],[1.0,0.0]], [1.0,1.0,1.0,0.01,0.0,0.0])");
  Register(
      loader, "chem_acquisition", Acquisition,
      "Expected improvement: rows [mean,sigma], parameters [best,minimize,xi].",
      "chem_acquisition([[2.0,1.0]], [2.0,0.0,0.0])");
  Register(loader, "chem_diversity_select", DiversitySelect,
           "Farthest-first candidate selection in supplied feature units: "
           "[count,existing_prefix_count].",
           "chem_diversity_select([[0.0],[1.0],[5.0],[10.0]], [2.0,1.0])");
  Register(loader, "chem_summarize", Summarize,
           "Column-wise descriptive statistics.",
           "chem_summarize([[1.0],[2.0],[3.0]], []::DOUBLE[])");
  Register(loader, "chem_welch_t_test", WelchTest,
           "Two-sided Welch t-test: rows [value,group 0/1].",
           "chem_welch_t_test([[1.0,0.0],[2.0,0.0],[3.0,0.0],[2.0,1.0],[3.0,1."
           "0],[4.0,1.0]], []::DOUBLE[])");
  Register(loader, "chem_one_way_anova", Anova,
           "One-way ANOVA: rows [value,group_id].",
           "chem_one_way_anova([[1.0,0.0],[2.0,0.0],[3.0,0.0],[2.0,1.0],[3.0,1."
           "0],[4.0,1.0]], []::DOUBLE[])");
  Register(loader, "chem_adjust_pvalues", AdjustPvalues,
           "Multiplicity adjustment: rows [p], parameters [0 Holm/1 BH].",
           "chem_adjust_pvalues([[0.01],[0.04],[0.03]], [1.0])");
  Register(loader, "chem_process_capability", ProcessCapability,
           "Overall sample-SD Pp/Ppk and descriptive mean +/-3SD limits: "
           "[LSL,USL,optional within_sigma]. Cp/Cpk require supplied "
           "within_sigma. Assumes an appropriate stable normal process.",
           "chem_process_capability([[1.0],[2.0],[3.0]], [-1.0,5.0])");
  Register(loader, "chem_kaplan_meier", KaplanMeier,
           "Kaplan-Meier survival for rows [time,event 0/1], with tied "
           "event/censor risk sets.",
           "chem_kaplan_meier([[1.0,1.0],[2.0,0.0],[3.0,1.0]], []::DOUBLE[])");
  Register(loader, "chem_random_intercept", RandomIntercept,
           "Balanced one-way random-intercept ANOVA variance estimator: rows "
           "[value,group_id]. Not a general mixed-effects solver.",
           "chem_random_intercept([[1.0,0.0],[3.0,0.0],[5.0,1.0],[7.0,1.0]], "
           "[]::DOUBLE[])");
}
} // namespace duckdb
