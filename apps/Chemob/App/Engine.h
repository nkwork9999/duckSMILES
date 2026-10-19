#include <stdint.h>
#include <stddef.h>
int32_t ds_workbench_analyze(const uint8_t *input, size_t length, uint8_t *output, size_t capacity);
int32_t ds_mol_graph_json(const uint8_t *input, size_t length, uint8_t *output, size_t capacity);
int32_t ds_workbench_drawing_graph(const uint8_t *input, size_t length, uint8_t *output, size_t capacity);
int32_t ds_workbench_revision(uint8_t *output, size_t capacity);

int32_t ds_workbench_catalog(uint8_t *output, size_t capacity);
