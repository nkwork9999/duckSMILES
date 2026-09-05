"use strict";

const assert = require("node:assert/strict");
const test = require("node:test");
const { checkWasmLinkage } = require("./check_wasm_linkage.cjs");

// Small real Wasm modules keep the validator tests independent of an SDK.
function moduleBytes(importName, entry = "ducksmiles_duckdb_cpp_init") {
  const text = (value) => [value.length, ...Buffer.from(value)];
  const section = (id, payload) => [id, payload.length, ...payload];
  return Buffer.from([
    0, 97, 115, 109, 1, 0, 0, 0,
    ...section(1, [1, 0x60, 0, 0]),
    ...(importName ? section(2, [1, ...text("env"), ...text(importName), 0, 0]) : []),
    ...section(3, [1, 0]),
    ...section(7, [1, ...text(entry), 0, importName ? 1 : 0]),
    ...section(10, [1, 2, 0, 0x0b]),
  ]);
}

test("accepts a linked module and ordinary host imports", () => {
  assert.doesNotThrow(() => checkWasmLinkage(moduleBytes()));
  assert.doesNotThrow(() => checkWasmLinkage(moduleBytes("malloc")));
});

test("rejects unresolved descriptor functions", () => {
  assert.throws(() => checkWasmLinkage(moduleBytes("ds_mol_num_atoms")), /Unlinked Rust/);
});

test("rejects unresolved Rust runtime symbols", () => {
  assert.throws(() => checkWasmLinkage(moduleBytes("_RNvCs_test")), /Unlinked Rust/);
  assert.throws(() => checkWasmLinkage(moduleBytes("_ZN4test17h0123456789abcdefE")), /Unlinked Rust/);
});

test("rejects missing entry point and malformed binaries", () => {
  assert.throws(() => checkWasmLinkage(moduleBytes(undefined, "other")), /Missing DuckSMILES/);
  assert.throws(() => checkWasmLinkage(Buffer.from("not a Wasm binary")), WebAssembly.CompileError);
});
