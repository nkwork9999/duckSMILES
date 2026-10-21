#pragma once
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
// Returns an owned UTF-8 JSON buffer, or NULL for invalid/numerically
// unsupported inputs. Only ds_analytics_free may release this exact
// pointer/length pair.
uint8_t *ds_analytics_run(const uint8_t *operation, size_t operation_len,
                          const double *data, size_t rows, size_t columns,
                          const double *parameters, size_t parameters_len,
                          size_t *output_len);
void ds_analytics_free(uint8_t *data, size_t len);
#ifdef __cplusplus
}
#endif
