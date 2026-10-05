# Improvement progress

Target: 100 independently reviewable improvements. Completed in this batch: 11. The target is not complete. Test cases and formatting edits are not counted separately.

1. XYZ atoms with invalid/non-finite coordinates are skipped rather than placed at zero.
2. Short or non-ASCII MODEL records no longer panic.
3. PDB atom record recognition excludes ATOMIC and other prefix collisions.
4. PDB atoms with invalid/non-finite coordinates are skipped.
5. mmCIF atom loops continue past blank/comment lines.
6. mmCIF parsing recognizes data/save/stop loop boundaries.
7. mmCIF rows with unsupported group_PDB values are excluded.
8. mmCIF atoms require finite x/y/z coordinates.
9. mmCIF missing metadata markers become empty fields.
10. mmCIF fields support both single and double quotes.
11. Unterminated mmCIF quoted fields are rejected.

## Validation

Locked workspace tests passed (including parser regressions; one existing ignored test). The affected pdb crate passes Clippy. Workspace Clippy fails in existing smiles FFI raw-pointer APIs; these APIs were not changed.
