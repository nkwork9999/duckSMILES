# Rebuild std with the same Emscripten ABI/features as the side module.
# Precompiled std contains compiler_builtins objects that reject shared memory.
# https://doc.rust-lang.org/rustc/platform-support/wasm32-unknown-emscripten.html
set(DUCKSMILES_WASM_RUST_TOOLCHAIN "nightly-2025-02-01" CACHE STRING
    "Pinned Rust toolchain used to rebuild the Wasm standard library")
find_program(DUCKSMILES_RUSTUP_EXECUTABLE rustup)
if(NOT DUCKSMILES_RUSTUP_EXECUTABLE)
    message(FATAL_ERROR "rustup is required to prepare the Wasm standard library")
endif()
execute_process(
    COMMAND "${DUCKSMILES_RUSTUP_EXECUTABLE}" toolchain install
            "${DUCKSMILES_WASM_RUST_TOOLCHAIN}" --profile minimal
            --component rust-src --no-self-update
    RESULT_VARIABLE DUCKSMILES_WASM_TOOLCHAIN_RESULT
)
if(NOT DUCKSMILES_WASM_TOOLCHAIN_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to prepare the pinned Wasm Rust toolchain")
endif()

set(CARGO_EXTRA_FLAGS -Zbuild-std=std,panic_abort
    -Zbuild-std-features=panic_immediate_abort)
set(CARGO_EXTRA_ENV
    "RUSTUP_TOOLCHAIN=${DUCKSMILES_WASM_RUST_TOOLCHAIN}"
    "CARGO_PROFILE_RELEASE_PANIC=abort"
    "CARGO_PROFILE_DEV_PANIC=abort"
)
if(USE_WASM_THREADS)
    list(APPEND CARGO_EXTRA_ENV
        "RUSTFLAGS=$ENV{RUSTFLAGS} -Ctarget-feature=+atomics,+bulk-memory,+mutable-globals"
        "EMCC_CFLAGS=$ENV{EMCC_CFLAGS} -pthread"
    )
endif()
