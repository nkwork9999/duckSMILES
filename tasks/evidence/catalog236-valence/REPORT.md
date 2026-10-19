# Partial aromatic valence and cyclic triple bonds — v25

Frozen host checkpoint only. 212/217 RDKit2026.03.6 standards +24profiles =236.
Five standards and additional50 remain unfinished at this checkpoint. The
installed iPad App remains v24; current v26 source is separate ongoing work.

## Correction and tests

Numeric aromatic bond valence, explicit-valence rounding and implicit H follow
pinned Atom.cpp. Prepared acyclic AROMATIC bond type need not have an aromatic
flag. Cyclic aromatic TRIPLE retains its type. Surface/TPSA, conjugation,
ring classes, weighted distance and Bertz consume the appropriate flag/type.
BCUT/Morgan retain bond-type semantics. Raw parser/native graph consumers retain
their existing grammar/H contract. No tolerance relaxed or runtime dependency.

Independent41 inputs:205 H/mass values+41formulas, plus all236 fields on18 small
partial-colon/triple controls. Runtime failures progressed85mass/H cells,
593expanded cells,203surface cells,21remaining cells, then60triple-ring cells.
All selected values pass; all41 per-atom H vectors pass. Full release workspace
772passed,0failed,1existingignored; zero warnings. Full C ABI413 controls pass.

All236 fields+formula on41 inputs:37fully match,4still differ (three porphyrins
and one mixed aromatic/double cycle). All41formulas match. These broader
aromaticity/Kekulization gaps are retained in v25-controls41, not hidden by the
narrow selected-field test. The selected test does not assert full41×236 parity.

## Wider difficult set

All415 prior v20 failed/unsupported inputs, no native acceptance filter:
21all236+formula match,146accepted rows have numeric differences,248unsupported.
Formula mismatches5→0, fully matching19→21; no introduced differing rows or
newly differing fields within this selected cohort. QED25→23. Strict gateFAIL.
This is not a new whole100k count or a claim of complete RDKit precision.

## Actual Swift I/O and cache

Actual App Chemistry.swift/Store on Mac, no rendered UI:10000inputs,25unsupported,
236numeric/240CSV columns,2352408numeric cells independently checked. Source
text/order/CSV reimport preserved, all236/base/formula/status exactly unchanged
fromv24 on this particular10k corpus. Ten thousand warm rows reused with zero
new native calculations. The held-out inputs, rather than this10k set, expose
these corrections. Full10k disk-capacity proof remainsv24; focused corrected
C:C/benzyne disk/cache harness failsv24 and passesv25. Revision namespaces differ.

Uncontrolled Mac observation: cold46.652s,warm.160s,export/snapshot13.704s.
Not a controlled version/RDKit comparison, physical-device or UI measurement.
Current iPad binary remainsv24 while the next50-descriptor source is developed.
