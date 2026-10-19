# v10 cyclic-donation experiment — superseded

Native revision `ducksmiles-workbench-v10-cyclic139-rdkit2026036`.
139 entries (115 standard + 24 profile), additional50 unimplemented.

- Rust 736 passed / zero failed / one ignored.
- NCI 4,999 rows: 4,991 both accepted, 4,986 all139+formula matching;
  five accepted value differences + eight acceptance differences. Strict failed.
- Complete ChEBI 100,000 distinct graphs: 99,755 common accepted,
  99,099 all139+formula matching; 656 accepted value differences + 245
  unsupported = 901 differing rows. Six mass/formula differences remain.
- Compared with v8, 324 rows resolved but **nine previously matching rows now
  differ**. `chebi/scope-delta.json` records exact source indices.
  This build was not reflected into the Simulator app; v11 adds candidate-only
  connected-subsystem detection. See `../catalog139-fused/REPORT.md`.
- `cache-bulk.json`: actual Swift ChemistryStore disk save/reopen at capacity,
  10,000 distinct inputs, 25 unsupported retained, 8,669,292-byte file;
  all 10,000 reused, zero native calls, all values and row positions equal.
- `unsupported-baseline-audit.json`: v8's 245 unsupported source rows consist
  of 102 above the atom bound, 107 CXSMILES and 36 non-tetrahedral annotations.
  Do not discard annotations and call them fully supported.

The strict gate failed. A lower aggregate error count did not hide the nine
new regressions. Reports and hashes are retained as an experiment, not a final
compatibility claim.
