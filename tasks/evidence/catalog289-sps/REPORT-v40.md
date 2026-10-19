# Native v40 checkpoint — expanded preparation domain

2026-10-18. Native revision `ducksmiles-workbench-v40-large291-rdkit2026036`.
Host dylib SHA `9a2607acd30540ab29ed41004e2cc118b03ab39ff28519d2e7dc82fc50bc23c5`.
The normal installed Simulator App is still v38. No v40 App GUI claim.

Prepared workbench accepts up to640atoms/1280bonds within4096input bytes. Raw
graph sidecar remains256/512. Six independent boundaries and the257-atom
all291+formula control pass after actual runtime RED. [Tests](large-input-tests.json).
Whole100k preparation now99995accepted/5unsupported,102new and zero regressions.
[All source positions](v40-acceptance100k.json). This is input support only.

All102 original size-limited molecules were compared using the frozen C ABI,
with no native selection. All accepted, all formulas agree;95 match all291.
Seven differ ONLY in Harary, absolute errors1.001e-9..1.439e-9 against the
unchanged1e-9 tolerance. No tolerance widening; checkpoint gate FAILS.
Independent exact rational distance sums confirm the dense-reference values
for these seven. A runtime regression test is being added before correction.
[Complete comparison](v40-large102/comparison.json),
[Exact rational references](large-next-controls/harary-exact-reference.json).
Pinned Ipc/AvgIpc compatibility differs on all102 separately; native stable
definitions agree with the independent exact coefficient reference.

Full release Rust799passed/0failed/2ignored/0warnings. An older integration
assertion still rejected257atoms; its expected boundary was updated to641,
with the failure preserved. [Workspace](v40-workspace.json).
An initial library-copy mistake loaded frozen39: its runtime revision exposed
the mistake, and that attempt is explicitly INVALID evidence under
`v40-stale-v39-library-rejected`. The correct static archive was then linked,
and runtime revision/hash validated before the complete102 comparison.

Actual Swift Store tests pass: real frozen38 negative geometry and257-atom
entries are recomputed, duplicate reuse2, disk reopen reuse4 with zero new
calculations and identical291-field results. All earlier cache regressions pass.
[Cache proof](v40-cache-integration.json). This is a host Store test, not GUI.

The598-atom maximum source molecule needs about11–16seconds on the loaded Mac
for all291 plus JSON, without cache or output retry. C-process peak RSS about
104MiB. Reference audit/Swift compilation overlapped, so this is a capacity
diagnostic, not a fair benchmark or physical iOS result. First3seconds of native
sampling are entirely in exact information coefficient computation.
[Measurement policy](v40-largest-measurement-policy.json).

Frozen37 remains the whole100k numerical baseline. Full latest-version100k
numerical rerun and final App integration are pending. Five multi-metal
canonical coordination cases remain unsupported. No completion claim.
