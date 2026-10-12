#pragma once

#include "duckdb/common/vector_operations/unary_executor.hpp"
#if __has_include("duckdb/common/vector/list_vector.hpp")
#include "duckdb/common/vector/list_vector.hpp"
#endif

namespace duckdb {
namespace ducksmiles_compat {

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

static void ToUnifiedFormat(Vector &input, idx_t count,
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

static auto ListChild(Vector &input) -> decltype(ListChildImpl(input, 0)) {
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
