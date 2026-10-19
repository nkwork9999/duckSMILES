"""Finalize descriptive catalog metadata after the frozen calculation audit.

Requires the original v42 source snapshot to match every calculation source.
Only CATALOG_JSON may change; the compiled SPECS prefix remains byte-identical.
Verification describes that calculation revision, never untested future changes.
"""
import argparse
import hashlib
import json
import re
import tarfile
from pathlib import Path


def practical_use(id):
    uses = {
        'ds_Wiener': 'Compare overall graph separation and branching in molecular feature tables, preferably within similar heavy-atom sizes; this is not a 3D distance.',
        'ds_Harary': 'Describe graph compactness through reciprocal topological distances for structural comparison; size effects and disconnected-component policy must be retained.',
        'ds_EccentricConnectivity': 'Describe how branching degree and graph reach combine at atoms, as a candidate structural feature for molecule cohorts.',
        'ds_DegreeDistance': 'Compare how high-degree sites are separated in the molecular graph; use alongside size and composition when comparing compounds.',
        'ds_Gutman': 'Describe degree-weighted graph separation for structural feature analysis, retaining the stated finite-component policy.',
    }
    if id in uses:
        return uses[id]
    m = re.fullmatch(r'ds_MDE([CNO])([1-4])([1-4])', id)
    if m:
        element = dict(C='carbon', N='nitrogen', O='oxygen')[m[1]]
        return (f'Compare separation of {element} sites with graph degrees {m[2]} and {m[3]} '
                'in structural feature tables. Degree is graph connectivity, not chemical valence; absent eligible pairs are undefined.')
    m = re.fullmatch(r'ds_GGI(\d+)', id)
    if m:
        return (f'Describe adjacency-distance asymmetry over atom pairs {m[1]} bonds apart as a topology feature. '
                'This is not an electrostatic or measured partial charge.')
    m = re.fullmatch(r'ds_MWC(\d+)', id)
    if m:
        return (f'Describe branching propagation through walks of length {m[1]} for molecular feature analysis. '
                'Vertices and edges may repeat; this is not a simple-path count.')
    m = re.fullmatch(r'ds_SRW(\d+)', id)
    if m:
        return (f'Describe returning connectivity through closed walks of length {m[1]}. '
                'Repeated edges are allowed, so even-order values can occur in trees; this is not a ring count.')
    raise AssertionError(id)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    report = json.loads(args.report.read_text())
    assert report['revision'] == 'ducksmiles-workbench-v42-cluster291-rdkit2026036'
    assert report['gate'] == report['accepted_value_gate'] == 'passed'
    assert report['all_values_matching_rows'] == report['common_valid'] == 100000
    catalog_path = root/'crates/smiles/src/descriptor_catalog.rs'
    marker = 'pub const CATALOG_JSON: &str = r###"'
    current = catalog_path.read_text()
    prefix, payload = current.split(marker, 1)
    raw, suffix = payload.split('"###;', 1)
    catalog = json.loads(raw)
    matched = {}
    ignored_sidecars = 0
    with tarfile.open(args.snapshot) as archive:
        for member in archive.getmembers():
            if not member.isfile() or not member.name.startswith('crates/smiles/src/') or not member.name.endswith('.rs'):
                continue
            if Path(member.name).name.startswith('._'):
                ignored_sidecars += 1
                continue  # macOS archive resource-fork metadata, not Rust source
            before = archive.extractfile(member).read()
            path = root/member.name
            if path == catalog_path:
                assert before.decode().split(marker,1)[0] == prefix
                baseline = json.loads(before.decode().split(marker,1)[1].split('"###;',1)[0])
                baseline_catalog_sha256 = hashlib.sha256(before).hexdigest()
            else:
                assert path.read_bytes() == before, member.name
                matched[member.name] = hashlib.sha256(before).hexdigest()
    assert len(matched) > 30 and len(catalog['descriptors']) == 291
    assert catalog == baseline, 'Start only from the frozen audit catalog'
    definitions_path = root/'tasks/descriptors/extra50-definitions.json'
    definitions = json.loads(definitions_path.read_text())
    extras = {row['id']: row for row in definitions['descriptors']}
    assert len(extras) == 50
    for row in catalog['descriptors']:
        if row['id'] in ('Ipc', 'AvgIpc'):
            row['verification'] = 'verified-exact-definition-v42-input100000; pinned numerical compatibility reported separately'
        else:
            row['verification'] = 'verified-declared-reference-v42-input100000'
        if row['origin'] == 'native-extra':
            row['practical_use'] = practical_use(row['id'])
            extra = extras[row['id']]
            extra['practical_use'] = row['practical_use']
            extra['verification'] = row['verification']
    def mathematical_rows(rows):
        return [{k:v for k,v in row.items() if k not in ('verification','practical_use')} for row in rows]
    assert mathematical_rows(catalog['descriptors']) == mathematical_rows(baseline['descriptors'])
    encoded = json.dumps(catalog, ensure_ascii=False, separators=(',', ':'))
    catalog_path.write_text(prefix + marker + encoded + '"###;' + suffix)
    assert catalog_path.read_text().split(marker,1)[0] == prefix
    definitions_path.write_text(json.dumps(definitions,ensure_ascii=False,indent=2)+'\n')
    phase_path = root/'tasks/descriptors/phase1-catalog.json'
    phase = json.loads(phase_path.read_text())
    assert mathematical_rows(phase['descriptors']) == mathematical_rows(baseline['descriptors'])
    phase['descriptors'] = catalog['descriptors']
    phase_path.write_text(json.dumps(phase,ensure_ascii=False,indent=2)+'\n')
    result = dict(calculation_revision=report['revision'],calculation_audit_library_sha256=report['library_sha256'],
                  audit_report_sha256=hashlib.sha256(args.report.read_bytes()).hexdigest(),
                  numerical_source_files_identical=matched,
                  ignored_archive_resource_fork_sidecars=ignored_sidecars,
                  baseline_catalog_source_sha256=baseline_catalog_sha256,
                  specs_prefix_sha256=hashlib.sha256(prefix.encode()).hexdigest(),
                  specs_prefix_byte_identical=True,all_catalog_math_order_ids_groups_units_digits_tolerances_unchanged=True,
                  changed_fields_only=['verification','practical_use'],practical_uses_defined=50,
                  new_catalog_source_sha256=hashlib.sha256(catalog_path.read_bytes()).hexdigest(),
                  calculation_namespace_change=False,predictive_performance_claim=False,
                  new_binary_build_pending=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='numerical_source_files_identical'}))


if __name__ == '__main__':
    main()
