#pragma once

#include "duckdb/common/vector_operations/unary_executor.hpp"
#if __has_include("duckdb/common/vector/list_vector.hpp")
#include "duckdb/common/vector/list_vector.hpp"
#endif

namespace duckdb {
namespace ducksmiles_compat {

// DuckDB 2 stores scalar argument types in its function signature. Preserve
// declaration order when associating catalog documentation with an overload.
template <class F>
static auto ArgumentTypesImpl(const F &function, int)
    -> decltype(function.GetSignature(), vector<LogicalType>()) {
  vector<LogicalType> types;
  for (const auto &parameter : function.GetSignature().GetParameters()) {
    types.push_back(parameter.GetType());
  }
  return types;
}

template <class F>
static auto ArgumentTypesImpl(const F &function, long)
    -> decltype(function.arguments) {
  return function.arguments;
}

template <class F> static vector<LogicalType> ArgumentTypes(const F &function) {
  return ArgumentTypesImpl(function, 0);
}

// DuckDB 2's catalog reads parameter names from the signature itself. Keep
// their positional-only kind so documenting names does not change binding.
template <class F>
static auto SetArgumentNamesImpl(F &function, const vector<string> &names, int)
    -> decltype(function.GetSignature(), void()) {
  for (idx_t i = 0; i < names.size(); i++) {
    function.GetSignature().GetParameter(i).SetName(names[i].c_str());
  }
}

template <class F>
static void SetArgumentNamesImpl(F &, const vector<string> &, long) {}

template <class F>
static void SetArgumentNames(F &function, const vector<string> &names) {
  SetArgumentNamesImpl(function, names, 0);
}

// DuckDB 2 keeps vector sizes internally. Select its count-free APIs without
// using deprecated overloads, while retaining the explicit count on DuckDB 1.
template <class V>
static auto FlattenImpl(V &input, idx_t, int)
    -> decltype(input.Flatten(), void()) {
  input.Flatten();
}

template <class V> static void FlattenImpl(V &input, idx_t count, long) {
  input.Flatten(count);
}

template <class V> static void Flatten(V &input, idx_t count) {
  FlattenImpl(input, count, 0);
}

template <class V>
static auto ToUnifiedFormatImpl(V &input, idx_t, UnifiedVectorFormat &data, int)
    -> decltype(input.ToUnifiedFormat(data), void()) {
  input.ToUnifiedFormat(data);
}

template <class V>
static void ToUnifiedFormatImpl(V &input, idx_t count,
                                UnifiedVectorFormat &data, long) {
  input.ToUnifiedFormat(count, data);
}

static inline void ToUnifiedFormat(Vector &input, idx_t count,
                            UnifiedVectorFormat &data) {
  ToUnifiedFormatImpl(input, count, data, 0);
}

template <class List = ListVector>
static auto ListChildImpl(Vector &input, int)
    -> decltype(List::GetChild(input)) {
  return List::GetChild(input);
}

template <class List = ListVector>
static Vector &ListChildImpl(Vector &input, long) {
  return List::GetEntry(input);
}

static inline auto ListChild(Vector &input) -> decltype(ListChildImpl(input, 0)) {
  return ListChildImpl(input, 0);
}

// GenericExecute is shared by stable DuckDB and the new scalar executor.
// Keep the callback's validity mask so Rust error sentinels still become NULL.
template <class Callback> struct NullableUnaryOperation {
  template <class Input, class Output>
  static Output Operation(Input input, ValidityMask &mask, idx_t index,
                          void *data) {
    return (*static_cast<Callback *>(data))(input, mask, index);
  }
};

template <class Input, class Output, class Callback>
static void ExecuteWithNulls(Vector &input, Vector &result, idx_t count,
                             Callback callback) {
  void *data = &callback;
  UnaryExecutor::GenericExecute<Input, Output,
                                NullableUnaryOperation<Callback>>(
      input, result, count, data, true);
}

// New DuckDB vectors expose read-only data by default. Use their mutable API
// when available, without casting away constness or modifying input buffers.
template <class T, class Flat = FlatVector>
static auto ResultDataImpl(Vector &result, idx_t count, int)
    -> decltype(Flat::template GetDataMutable<T>(result)) {
  Flat::SetSize(result, count);
  return Flat::template GetDataMutable<T>(result);
}

template <class T, class Flat = FlatVector>
static T *ResultDataImpl(Vector &result, idx_t count, long) {
  return Flat::template GetData<T>(result);
}

template <class T> static T *ResultData(Vector &result, idx_t count) {
  return ResultDataImpl<T>(result, count, 0);
}

template <class Flat = FlatVector>
static auto ResultValidityImpl(Vector &result, int)
    -> decltype(Flat::ValidityMutable(result)) {
  return Flat::ValidityMutable(result);
}

template <class Flat = FlatVector>
static ValidityMask &ResultValidityImpl(Vector &result, long) {
  return Flat::Validity(result);
}

static ValidityMask &ResultValidity(Vector &result) {
  return ResultValidityImpl(result, 0);
}

} // namespace ducksmiles_compat
} // namespace duckdb
