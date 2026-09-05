"use strict";

const fs = require("node:fs");

function checkWasmLinkage(bytes) {
  const module = new WebAssembly.Module(bytes);
  const exports = WebAssembly.Module.exports(module);
  const unresolved = WebAssembly.Module.imports(module).filter((symbol) => {
    if (!/^ds_|^_R|^_ZN.*17h[0-9a-f]{16}E$|^__(rust|rdl|rg)_|^rust_(eh_personality|begin_unwind)$/.test(symbol.name)) {
      return false;
    }
    // Emscripten's loader fills GOT entries from the side module's own exports.
    // These relocations are not missing Rust definitions. Direct env imports,
    // missing exports and mismatched export kinds must still fail.
    const exportKind = symbol.module === "GOT.func" ? "function" :
      symbol.module === "GOT.mem" ? "global" : undefined;
    return !(symbol.kind === "global" && exportKind && exports.some(
      ({ name, kind }) => name === symbol.name && kind === exportKind
    ));
  });
  if (unresolved.length) {
    throw new Error(
      `Unlinked Rust symbols (${unresolved.length}): ` +
        unresolved.slice(0, 10).map(({ name }) => name).join(", ")
    );
  }
  const hasEntry = exports.some(
    ({ name, kind }) => name === "ducksmiles_duckdb_cpp_init" && kind === "function"
  );
  if (!hasEntry) {
    throw new Error("Missing DuckSMILES extension entry point");
  }
}

module.exports = { checkWasmLinkage };

if (require.main === module) {
  if (process.argv.length !== 3) {
    throw new Error("Usage: node scripts/check_wasm_linkage.cjs <extension.wasm>");
  }
  checkWasmLinkage(fs.readFileSync(process.argv[2]));
  console.log("Wasm linkage verified: entry point present; Rust imports resolved within side module");
}
