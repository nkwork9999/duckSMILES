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

function gotModuleBytes(namespace, exportedKind) {
  const name = "__rust_alloc";
  const text = (value) => [value.length, ...Buffer.from(value)];
  const section = (id, payload) => [id, payload.length, ...payload];
  return Buffer.from([
    0, 97, 115, 109, 1, 0, 0, 0,
    ...section(1, [1, 0x60, 0, 0]),
    ...section(2, [1, ...text(namespace), ...text(name), 3, 0x7f, 1]),
    ...section(3, [1, 0]),
    ...section(6, [1, 0x7f, 0, 0x41, 1, 0x0b]),
    ...section(7, [exportedKind === undefined ? 1 : 2,
      ...text("ducksmiles_duckdb_cpp_init"), 0, 0,
      ...(exportedKind === undefined ? [] : [
        ...text(name), exportedKind, exportedKind === 3 ? 1 : 0,
      ]),
    ]),
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
  assert.throws(() => checkWasmLinkage(moduleBytes("__rust_alloc")), /Unlinked Rust/);
  assert.throws(() => checkWasmLinkage(moduleBytes("rust_eh_personality")), /Unlinked Rust/);
});

test("accepts GOT relocations resolved by the side module's own exports", () => {
  assert.doesNotThrow(() => checkWasmLinkage(gotModuleBytes("GOT.func", 0)));
  assert.doesNotThrow(() => checkWasmLinkage(gotModuleBytes("GOT.mem", 3)));
});

test("rejects missing or mismatched GOT definitions and direct env imports", () => {
  for (const bytes of [
    gotModuleBytes("GOT.func"), gotModuleBytes("GOT.mem"),
    gotModuleBytes("GOT.func", 3), gotModuleBytes("GOT.mem", 0),
    gotModuleBytes("env", 3),
  ]) {
    assert.throws(() => checkWasmLinkage(bytes), /Unlinked Rust/);
  }
});

test("rejects missing entry point and malformed binaries", () => {
  assert.throws(() => checkWasmLinkage(moduleBytes(undefined, "other")), /Missing DuckSMILES/);
  assert.throws(() => checkWasmLinkage(Buffer.from("not a Wasm binary")), WebAssembly.CompileError);
});
