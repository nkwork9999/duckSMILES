// macOS C ABI capacity measurement. No Python interpreter in the measured process.
// Usage: measure_one_native library input.smi result.json [calls1..4]
// One cold and three repeated computations; no result cache or buffer retry.
#include <assert.h>
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>

typedef int32_t (*analyze_fn)(const uint8_t *, size_t, uint8_t *, size_t);
typedef int32_t (*revision_fn)(uint8_t *, size_t);

static double now(void) {
  struct timespec value;
  assert(clock_gettime(CLOCK_MONOTONIC, &value) == 0);
  return (double)value.tv_sec + (double)value.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
  assert(argc == 4 || argc == 5);
  int calls = argc == 5 ? atoi(argv[4]) : 4;
  assert(calls >= 1 && calls <= 4);
  void *library = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
  if (!library) { fprintf(stderr, "%s\n", dlerror()); return 2; }
  analyze_fn analyze = (analyze_fn)dlsym(library, "ds_workbench_analyze");
  revision_fn revision = (revision_fn)dlsym(library, "ds_workbench_revision");
  assert(analyze && revision);
  uint8_t version[256] = {0};
  int32_t vlen = revision(version, sizeof version - 1);
  assert(vlen > 0 && vlen < (int32_t)sizeof version);
  assert(strcmp((char *)version, "ducksmiles-workbench-v41-harary291-rdkit2026036") == 0);
  FILE *source = fopen(argv[2], "rb");
  assert(source);
  uint8_t input[4097];
  size_t length = fread(input, 1, sizeof input, source);
  assert(!ferror(source) && fclose(source) == 0);
  while (length && (input[length - 1] == '\n' || input[length - 1] == '\r')) --length;
  assert(length && length <= 4096);
  const size_t capacity = 1024 * 1024;
  uint8_t *out = malloc(capacity);
  assert(out);
  double seconds[4];
  int32_t required = 0;
  for (int i = 0; i < calls; ++i) {
    double start = now();
    required = analyze(input, length, out, capacity);
    seconds[i] = now() - start;
    assert(required > 0 && required <= (int32_t)capacity);
  }
  FILE *result = fopen(argv[3], "wb");
  assert(result && fwrite(out, 1, (size_t)required, result) == (size_t)required);
  assert(fclose(result) == 0);
  struct rusage usage;
  assert(getrusage(RUSAGE_SELF, &usage) == 0);
  printf("{\"revision\":\"%s\",\"scope\":\"macOS C ABI, all291, no cache\","
         "\"cold_seconds\":%.9f,\"repeat_seconds\":[",
         version, seconds[0]);
  for (int i = 1; i < calls; ++i) printf("%s%.9f", i == 1 ? "" : ",", seconds[i]);
  printf("],"
         "\"process_max_rss_bytes\":%ld,\"output_bytes\":%d,\"input_bytes\":%zu,"
         "\"physical_ios_measurement\":false}\n",
         usage.ru_maxrss, required, length);
  free(out);
  assert(dlclose(library) == 0);
  return 0;
}
