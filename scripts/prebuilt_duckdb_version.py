#!/usr/bin/env python3
"""Align native prebuilt metadata without disabling DuckDB's ABI/version checks.

The official same-commit staging record disambiguates runtime versions from
storage-version history in the archive. Embedded strings are checked, rather
than interpreting the LibraryVersion function in ELF/COFF/Mach-O. Actual SQL
loading remains the final compatibility check.
"""

import argparse
import mmap
from pathlib import Path
import re
import subprocess
import sys
import urllib.request


VERSION = re.compile(r"v([0-9]+)\.([0-9]+)\.([0-9]+)(?:-(alpha|rc|dev)([0-9]+))?")
EMBEDDED_VERSION = re.compile(
    rb"(?:(?<=\x00)|\A)(v[0-9]+\.[0-9]+\.[0-9]+(?:-(?:alpha|rc|dev)[0-9]+)?)\x00"
)
MAX_LIBRARY_BYTES = 1024 * 1024 * 1024


def _base(value, description):
    match = VERSION.fullmatch(value)
    if not match:
        raise ValueError("invalid " + description + " version")
    return match.group(1, 2, 3)


def _contains_string(library, value):
    token = value.encode("ascii") + b"\0"
    start = 0
    while True:
        position = library.find(token, start)
        if position < 0:
            return False
        if position == 0 or library[position - 1] == 0:
            return True
        start = position + 1


def _check_source(source_id, library):
    if not re.fullmatch(r"[0-9a-f]{10}", source_id):
        raise ValueError("invalid source checkout Source ID")
    if not _contains_string(library, source_id):
        raise ValueError("Source ID does not match the prebuilt library")


def _versions(library):
    return {match.group(1).decode("ascii") for match in EMBEDDED_VERSION.finditer(library)}


def verify_version(source_id, requested_version, official_version, library):
    """Verify a record against literal binary strings and the requested base."""
    requested_base = _base(requested_version, "requested")
    _check_source(source_id, library)
    candidates = {value for value in _versions(library) if _base(value, "embedded") == requested_base}
    if official_version is None:
        if any("-" in value for value in candidates):
            raise ValueError("an official staging record is required for prerelease metadata")
        if requested_version not in candidates:
            raise ValueError("requested version is not embedded in the prebuilt library")
        return requested_version
    official_version = official_version.strip()
    official_base = _base(official_version, "official staging record")
    if official_base != requested_base:
        raise ValueError("official staging record has a different requested base version")
    if official_version not in candidates:
        raise ValueError("official staging record version is not embedded in the prebuilt library")
    return official_version


def _source_id(source_dir):
    source_dir = Path(source_dir).resolve()
    try:
        top = subprocess.check_output(
            ["git", "-C", str(source_dir), "rev-parse", "--show-toplevel"],
            stderr=subprocess.PIPE,
            text=True,
        ).strip()
        if Path(top).resolve() != source_dir:
            raise ValueError("source directory must be the DuckDB checkout root")
        commit = subprocess.check_output(
            ["git", "-C", str(source_dir), "rev-parse", "--verify", "HEAD"],
            stderr=subprocess.PIPE,
            text=True,
        ).strip().lower()
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError("cannot read the DuckDB source checkout HEAD") from error
    if not re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}", commit):
        raise ValueError("invalid DuckDB source checkout HEAD")
    return commit[:10]


def _fetch_record(source_id):
    request = urllib.request.Request(
        "https://duckdb-staging.duckdb.org/" + source_id + "/latest_alpha_version.txt",
        headers={"User-Agent": "duckdb-extension-ci-tools/prebuilt-duckdb"},
    )
    try:
        with urllib.request.urlopen(request, timeout=15) as response:
            data = response.read(4097)
        if len(data) > 4096:
            raise ValueError("official staging record is too large")
        return data.decode("ascii")
    except (OSError, UnicodeError) as error:
        raise ValueError("cannot fetch the official same-source staging record") from error


def resolve_version(library, source_dir, requested_version):
    source_id = _source_id(source_dir)
    requested_base = _base(requested_version, "requested")
    with Path(library).open("rb") as stream:
        size = stream.seek(0, 2)
        if not 0 < size <= MAX_LIBRARY_BYTES:
            raise ValueError("prebuilt library must contain 1 byte to 1 GiB")
        with mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_READ) as data:
            _check_source(source_id, data)
            prerelease = any(
                "-" in value and _base(value, "embedded") == requested_base
                for value in _versions(data)
            )
            official = _fetch_record(source_id) if prerelease else None
            return verify_version(source_id, requested_version, official, data)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", required=True, type=Path)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--requested-version", required=True)
    args = parser.parse_args(argv)
    try:
        result = resolve_version(args.library, args.source_dir, args.requested_version)
    except (OSError, ValueError) as error:
        print("prebuilt DuckDB version: " + str(error), file=sys.stderr)
        return 1
    print(result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
