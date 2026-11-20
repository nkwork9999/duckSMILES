"""Compare every profile coordinate in the Lean, Rust and optional DuckDB paths."""
import argparse
import csv
import json
import subprocess
import tempfile
import time
from pathlib import Path

from reference import CATALOG, CORPUS, observations, oracle

HERE = Path(__file__).resolve().parent


def run(library, extension, output, coverage):
    native = oracle.Native(library.resolve())
    specs = json.loads(CATALOG.read_text())["descriptors"]
    if native.catalog != specs:
        raise ValueError("native catalog or tolerances differ from reference")
    with tempfile.TemporaryDirectory() as tmp:
        ref, actual = (Path(tmp) / n for n in ("observations.json", "actual.json"))
        obs = observations(CORPUS)
        ref.write_text(json.dumps(obs, allow_nan=False))
        executable = HERE / ".lake/build/bin/descriptor_check"
        def lean(mode, *args, **kwargs):
            return subprocess.run([executable, mode, CATALOG, ref, *args],
                                  text=True, **kwargs)
        def check(rows):
            actual.write_text(json.dumps(rows, allow_nan=False))
            lean("check", actual, check=True)
        guide = [{"smiles": s, "values": [o["oracle"][c["id"]] for c in specs] if o else None}
                 for s, o in ((row["smiles"], row["observations"]) for row in obs)]
        check(guide)
        def rust():
            return [{"smiles": s, "values": row.get("descriptor_values")}
                    for s, row in zip(CORPUS, native.analyze(CORPUS), strict=True)]
        check(rust())
        if extension:
            import duckdb
            db = duckdb.connect(config={"allow_unsigned_extensions": "true"})
            path = str(extension.resolve()).replace("'", "''")
            db.execute(f"LOAD '{path}'")
            assert json.loads(db.execute("SELECT mol_descriptor_catalog_json()").fetchone()[0])["descriptors"] == specs
            def sql():
                raw = db.execute("SELECT s, mol_descriptor_values_json(s) FROM unnest(?::VARCHAR[]) WITH ORDINALITY t(s, i) ORDER BY i", [CORPUS]).fetchall()
                return [{"smiles": s, "values": json.loads(v) if v is not None else None} for s, v in raw]
            check(sql())
        if coverage:
            calculated = set(native.catalog[i]["id"] for i in range(len(specs))
                             if specs[i]["origin"] == "native-profile")
            calculated.update(("HeavyAtomCount", "NumHeteroatoms", "NOCount", "RingCount"))
            with coverage.open("w", newline="") as f:
                writer = csv.writer(f, lineterminator="\n")
                writer.writerow(["id", "reference", "Lean evaluation", "absolute tolerance", "relative tolerance",
                                 "defined reference values", "undefined reference values", "paths checked"])
                supported = [row["values"] for row in guide if row["values"] is not None]
                for i, spec in enumerate(specs):
                    key = spec["id"]
                    evaluator = ("atom/bond/ring aggregation" if key in calculated else
                                 "count unique SMARTS embeddings" if key.startswith("fr_") else
                                 "read declared numerical oracle")
                    defined = sum(row[i] is not None for row in supported)
                    writer.writerow([key, oracle.reference_source(spec), evaluator, spec["tolerance"],
                                     spec["relative_tolerance"], defined, len(supported) - defined,
                                     "Lean + Rust" + (" + SQL" if extension else "")])
        if output:
            cases = {"rust_ffi_batch": rust,
                     "lean_prepared_reference_process": lambda: json.loads(lean("reference", check=True, capture_output=True).stdout)}
            if extension:
                cases["duckdb_sql_batch"] = sql
            with output.open("w", newline="") as f:
                writer = csv.writer(f, lineterminator="\n")
                writer.writerow(["trial", "position", "path", "nanoseconds", "inputs"])
                for trial in range(6):
                    order = list(cases)
                    order = order[trial % len(order):] + order[:trial % len(order)]
                    for position, name in enumerate(order):
                        start = time.perf_counter_ns()
                        rows = cases[name]()
                        elapsed = time.perf_counter_ns() - start
                        check(rows)
                        writer.writerow([trial + 1, position + 1, name, elapsed, len(rows)])
        # Each mutation must fail, independently of every other mutation.
        for mutation in ("missing-column", "wrong-value", "undefined", "reorder-input", "lost-row"):
            rows = rust()
            if mutation == "missing-column":
                rows[0]["values"].pop()
            elif mutation == "wrong-value":
                rows[0]["values"][0] += 1
            elif mutation == "undefined":
                rows[0]["values"][6] = None
            elif mutation == "reorder-input":
                rows.reverse()
            else:
                rows.pop()
            actual.write_text(json.dumps(rows))
            res = lean("check", actual, capture_output=True)
            if res.returncode == 0:
                raise AssertionError(f"accepted {mutation}")
        print("Rejected all five independently mutated profiles.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("library", type=Path)
    parser.add_argument("--extension", type=Path)
    parser.add_argument("--timings", type=Path)
    parser.add_argument("--coverage", type=Path)
    args = parser.parse_args()
    run(args.library, args.extension, args.timings, args.coverage)
