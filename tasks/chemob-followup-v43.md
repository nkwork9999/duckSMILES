# Chemob v43 — checkpoint 2026-10-18

Local verification finished; physical iPhone execution and speed require a connected device.

- Prepared drawing API implemented; chemistry10,003 rows/RDKit0differences;14 size rejections.
- Original10k drawing coverage9,374;612 quality rejections;14 size rejections.261 partial stereo notices.
- Rust811 and Swift Core33 tests pass. Independent TH39,089 / E/Z3,593 round trip0differences.
- Simulator and unsigned iPhone arm64 Release builds pass.
- Compact iPhone chemical/stereo/column UI and iPad portrait/landscape inline UI pass.
- Original Japanese iPad six saved inputs retained; normal launch shows six reused results.
- Initial compact harness failure retained; final correction scrolls before reading offscreen rows.
- No descriptor arithmetic changes. Ipc/AvgIpc5723 numerical differences remain explicit stable columns.
- Physical phones unavailable; no physical speed claim. No commit, push or release.

Report: `tasks/evidence/depiction-ipc-v43/REPORT-v43.md`.
Machine checkpoint: `tasks/evidence/depiction-ipc-v43/resume-checkpoint.json`.
