"use strict";

const fs = require("node:fs");

function checkWasmLinkage(bytes) {
  const module = new WebAssembly.Module(bytes);
  const unresolved = WebAssembly.Module.imports(module).filter(({ name }) =>
    /^ds_|^_R|^_ZN.*17h[0-9a-f]{16}E$/.test(name)
  );
  if (unresolved.length) {
    throw new Error(
      `Unlinked Rust functions (${unresolved.length}): ` +
        unresolved.slice(0, 10).map(({ name }) => name).join(", ")
    );
  }
  const hasEntry = WebAssembly.Module.exports(module).some(
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
  console.log("Wasm linkage verified: extension entry point present; no Rust function imports");
}
