#pragma once

#include "duckdb/common/vector_operations/unary_executor.hpp"

namespace duckdb {
namespace ducksmiles_compat {

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
