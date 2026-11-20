"""Alternate the compiled programs on the same labels. Parsing precedes both timers."""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path
from reference import CORPUS, oracle

here = Path(__file__).resolve().parent
with tempfile.TemporaryDirectory() as tmp:
    atoms = Path(tmp) / "atoms"
    rows = []
    for smiles in CORPUS:
        m = oracle.Chem.MolFromSmiles(smiles, oracle.PARAMS) if smiles else None
        if m:
            rows.append([a.GetSymbol() for a in m.GetAtoms()])
    atoms.write_text("".join(" ".join(row) + "\n" for row in rows))
    rust = Path(tmp) / "counts"
    subprocess.run(["rustc", "-O", here / "bench_counts.rs", "-o", rust], check=True)
    programs = [("Lean", here / ".lake/build/bin/bench_counts"), ("Rust", rust)]
    repetitions = 100_000
    expected = repetitions // 5 * sum(a in ("C", "N", "O", "Cl", "F") for row in rows for a in row)
    out = csv.writer(sys.stdout, lineterminator="\n")
    out.writerow(["trial", "position", "language", "nanoseconds", "checksum"])
    for trial in range(6):
        for i, (name, executable) in enumerate(programs[::1 if trial % 2 == 0 else -1]):
            result = subprocess.check_output([executable, atoms, str(repetitions)], text=True)
            elapsed, checksum = map(int, result.strip().split(","))
            assert checksum == expected, (name, checksum, expected)
            out.writerow([trial + 1, i + 1, name, elapsed, checksum])
