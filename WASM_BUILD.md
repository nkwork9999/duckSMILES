# WebAssembly build contract

DuckDB builds Wasm extensions as side modules. All five Rust static archives
must be passed to DuckDB's final `emcc` command through
`DUCKDB_EXTENSION_DUCKSMILES_LINKED_LIBS`; `target_link_libraries` alone only
updates the intermediate C++ archive target.

The CI Emscripten SDK is **3.1.71**. Wasm builds use **nightly-2025-02-01**
(Rust 1.86 nightly / LLVM 19), with `rust-src`, to rebuild the Rust standard
library using `-Zbuild-std=std,panic_abort` and
`-Zbuild-std-features=panic_immediate_abort`. The Wasm-only Cargo profile uses
`panic=abort`, so Rust panics do not unwind across the C ABI. Native builds keep
their existing toolchain and panic behavior. No chemistry runtime dependency
is added.

For `wasm_threads`, Rust and its standard library are compiled with
`+atomics,+bulk-memory,+mutable-globals`, and Emscripten with `-pthread`.
Precompiled Rust standard libraries can be incompatible with these features
or the SDK's ABI; see the [Rust Emscripten target documentation](https://doc.rust-lang.org/rustc/platform-support/wasm32-unknown-emscripten.html).
Emscripten still labels side modules with pthreads experimental. A successful
build does not remove that upstream limitation.

`cmake/ConfigureRustWasm.cmake` installs the pinned build toolchain via rustup.
It can be selected explicitly with `DUCKSMILES_WASM_RUST_TOOLCHAIN`, but changing
that pin requires verifying all three Wasm variants against the selected SDK.
Use a clean CMake build directory and Cargo target directory when changing
SDK, Rust compiler, or thread/exception settings.

Every Wasm build runs `scripts/check_wasm_linkage.cjs` on the actual binary.
It checks the extension entry point and rejects unresolved DuckSMILES/Rust
symbols. Emscripten GOT relocations are accepted only when the side module
exports a matching function or global; its loader resolves those entries from
the module itself. The validator tests cover missing and mismatched definitions,
ordinary host imports, malformed binaries, and missing entry points:

```sh
node --test scripts/check_wasm_linkage.test.cjs
node scripts/check_wasm_linkage.cjs path/to/ducksmiles.duckdb_extension.wasm
```

This is a link-integrity check, not a DuckDB-Wasm runtime SQL test or a claim
of full RDKit compatibility. Molecular parity scope is documented separately
in [PROFILE_PARITY.md](PROFILE_PARITY.md).
