"""Compare every profile coordinate in the Lean, Rust and optional DuckDB paths."""
import argparse
import csv
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from reference import CATALOG, CORPUS, observations, oracle, oracle_values

HERE = Path(__file__).resolve().parent


def coverage_details(spec):
    k = spec["id"]
    finite = "catalog comparison only"
    if k.startswith("fr_") or k in ("NumHDonors", "NumHAcceptors", "NumAmideBonds", "NumRotatableBonds"):
        return "SMARTS search on graph and perceived atom properties", "fragmentCount_spec; mem_embeddings_iff; rooted_iff"
    if k.startswith(("SMR_VSA", "SlogP_VSA", "EState_VSA", "VSA_EState", "PEOE_VSA")):
        if k.startswith(("EState_VSA", "VSA_EState")):
            return "interval sums using calculated EState and calculated Labute areas", "accepts_shortest; bin_at_first_larger"
        return "interval sums using calculated Labute areas and prepared atomic contributions", "bin_at_first_larger"
    if k.startswith(("ds_MWC", "ds_SRW")):
        return "integer walk multiplicities on heavy graph; logarithmic reduction", "allWalks_counts; closedWalks_counts"
    if spec["origin"] == "native-extra":
        return "certified heavy graph distances and pair reductions", "accepts_shortest; mem_finitePairs"
    if spec["origin"] == "native-profile":
        return "integer reductions of perceived atom, bond and ring properties", "values_reordered; maximum_le"
    if k in ("NumHeterocycles",) or k.startswith(("NumAromatic", "NumAliphatic", "NumSaturated")):
        return "classify perceived ring basis", "classCounts_partition"
    if k in ("Ipc", "AvgIpc"):
        return "exact integer recurrence and binary64 entropy", "divideExact_correct (trace division)"
    if k.startswith("BCUT2D_"):
        return "Burden matrix and Jacobi eigensolver from hydrogen-depleted atom properties", finite
    if k == "qed":
        return "eight components computed from prepared graph and atom contributions; ADS and weighted geometric mean", "fragmentCount_spec; hasEmbedding_iff; accepts_shortest"
    if k.startswith("FpDensityMorgan"):
        return "initial connectivity hash and three rounds of Morgan refinement", finite
    if k.startswith(("Chi", "Kappa")) or k in ("BalabanJ", "BertzCT", "Phi"):
        return "graph connectivity, path or distance formula", finite
    if k == "LabuteASA":
        return "bond overlaps and periodic-table radii", finite
    if k.endswith("EStateIndex"):
        return "intrinsic states and distance perturbations", "accepts_shortest"
    if k in ("TPSA", "MolLogP", "MolMR", "HallKierAlpha") or k.endswith("PartialCharge"):
        return "reduction of prepared atomic contributions", finite
    return "atom, bond, stereo or perceived-ring aggregation", finite


def run(library, extension, output, coverage, extended=False):
    inputs = CORPUS
    if extended:
        sys.path.insert(0, str(HERE.parents[1] / "scripts"))
        from gen_rdkit_fixtures import MOLECULES
        inputs = list(dict.fromkeys(MOLECULES + ["[H]", "[H][H]", "[2H]O[2H]", "[NH4+]", "[Xe]", "C1", ""]))
    native = oracle.Native(library.resolve())
    specs = json.loads(CATALOG.read_text())["descriptors"]
    if native.catalog != specs:
        raise ValueError("native catalog or tolerances differ from reference")
    with tempfile.TemporaryDirectory() as tmp:
        ref, actual = (Path(tmp) / n for n in ("observations.json", "actual.json"))
        obs = observations(inputs)
        ref.write_text(json.dumps(obs, allow_nan=False))
        executable = HERE / ".lake/build/bin/descriptor_check"
        def lean(mode, *args, **kwargs):
            return subprocess.run([executable, mode, CATALOG, ref, *args],
                                  text=True, **kwargs)
        def check(rows):
            actual.write_text(json.dumps(rows, allow_nan=False))
            lean("check", actual, check=True)
        guide = [{"smiles": s, "values": oracle_values(s, specs)} for s in inputs]
        check(guide)
        def rust():
            return [{"smiles": s, "values": row.get("descriptor_values")}
                    for s, row in zip(inputs, native.analyze(inputs), strict=True)]
        check(rust())
        if extension:
            import duckdb
            db = duckdb.connect(config={"allow_unsigned_extensions": "true"})
            path = str(extension.resolve()).replace("'", "''")
            db.execute(f"LOAD '{path}'")
            assert json.loads(db.execute("SELECT mol_descriptor_catalog_json()").fetchone()[0])["descriptors"] == specs
            def sql():
                raw = db.execute("SELECT s, mol_descriptor_values_json(s) FROM unnest(?::VARCHAR[]) WITH ORDINALITY t(s, i) ORDER BY i", [inputs]).fetchall()
                return [{"smiles": s, "values": json.loads(v) if v is not None else None} for s, v in raw]
            check(sql())
        if coverage:
            with coverage.open("w", newline="") as f:
                writer = csv.writer(f, lineterminator="\n")
                writer.writerow(["id", "reference", "Lean evaluation", "general kernel result used", "absolute tolerance", "relative tolerance",
                                 "defined reference values", "undefined reference values", "paths checked"])
                supported = [row["values"] for row in guide if row["values"] is not None]
                for i, spec in enumerate(specs):
                    key = spec["id"]
                    evaluator, kernel = coverage_details(spec)
                    defined = sum(row[i] is not None for row in supported)
                    writer.writerow([key, oracle.reference_source(spec), evaluator, kernel, spec["tolerance"],
                                     spec["relative_tolerance"], defined, len(supported) - defined,
                                     "Lean + Rust" + (" + SQL" if extension else "")])
        if extended:
            inputs = CORPUS
            ref.write_text(json.dumps(observations(inputs), allow_nan=False))
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
    parser.add_argument("--extended", action="store_true", help="check the complete hand-curated molecule corpus")
    args = parser.parse_args()
    run(args.library, args.extension, args.timings, args.coverage, args.extended)
