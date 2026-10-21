#!/usr/bin/env python3
"""Fail-closed metadata checks for the native prebuilt build boundary."""

import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import prebuilt_duckdb_version as version


SOURCE_ID = "d54cbbe3df"


def library_bytes(*versions, source_id=SOURCE_ID):
    return b"archive\0" + b"\0".join(
        value.encode("ascii") for value in (*versions, source_id)
    ) + b"\0"


class VersionVerificationTests(unittest.TestCase):
    def test_official_record_selects_runtime_among_storage_history(self):
        data = library_bytes("v2.0.0", "v2.0.0-alpha41", "v2.0.0-alpha42")
        self.assertEqual(
            version.verify_version(SOURCE_ID, "v2.0.0", "v2.0.0-alpha42\n", data),
            "v2.0.0-alpha42",
        )

    def test_stable_only_library_needs_no_staging_record(self):
        data = library_bytes("v1.5.5", "v2.0.0", "v3.0.0-alpha1")
        self.assertEqual(version.verify_version(SOURCE_ID, "v2.0.0", None, data), "v2.0.0")

    def test_prerelease_without_official_record_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "official.*record"):
            version.verify_version(SOURCE_ID, "v2.0.0", None, library_bytes("v2.0.0", "v2.0.0-alpha42"))

    def test_record_must_have_same_requested_base(self):
        with self.assertRaisesRegex(ValueError, "base"):
            version.verify_version(SOURCE_ID, "v2.0.0", "v3.0.0-alpha1", library_bytes("v3.0.0-alpha1"))

    def test_record_must_exist_as_complete_embedded_string(self):
        with self.assertRaisesRegex(ValueError, "embedded"):
            version.verify_version(SOURCE_ID, "v2.0.0", "v2.0.0-alpha43", library_bytes("v2.0.0-alpha42"))
        data = b"prefixv2.0.0-alpha42\0" + SOURCE_ID.encode() + b"\0"
        with self.assertRaisesRegex(ValueError, "embedded"):
            version.verify_version(SOURCE_ID, "v2.0.0", "v2.0.0-alpha42", data)

    def test_source_id_must_match_complete_embedded_string(self):
        with self.assertRaisesRegex(ValueError, "Source ID"):
            version.verify_version(SOURCE_ID, "v2.0.0", "v2.0.0-alpha42", library_bytes("v2.0.0-alpha42", source_id="abcdef0123"))
        with self.assertRaisesRegex(ValueError, "Source ID"):
            version.verify_version(SOURCE_ID, "v2.0.0", None, b"v2.0.0\0prefix" + SOURCE_ID.encode() + b"\0")

    def test_invalid_or_multiple_record_values_are_rejected(self):
        for official in ["v2.0.0-alpha42\nv2.0.0-alpha43", "not-a-version", "v2.0.0-alpha42; echo unsafe"]:
            with self.subTest(official=official), self.assertRaisesRegex(ValueError, "record"):
                version.verify_version(SOURCE_ID, "v2.0.0", official, library_bytes("v2.0.0-alpha42"))

    def test_missing_version_and_invalid_request_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "version"):
            version.verify_version(SOURCE_ID, "v2.0.0", None, library_bytes("v1.5.5"))
        with self.assertRaisesRegex(ValueError, "requested"):
            version.verify_version(SOURCE_ID, "v2.0.0; echo unsafe", None, library_bytes("v2.0.0"))


class BuildBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="ducksmiles-prebuilt-test-")
        self.root = Path(self.temp.name)
        self.source = self.root / "source checkout"
        self.source.mkdir()
        subprocess.run(["git", "init", "-q", str(self.source)], check=True)
        subprocess.run(["git", "-C", str(self.source), "-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "commit", "-qm", "fixture", "--allow-empty"], check=True)
        self.source_id = subprocess.check_output(["git", "-C", str(self.source), "rev-parse", "HEAD"], text=True).strip()[:10]
        self.library = self.root / "lib 'quoted $cash `literal`.a"

    def tearDown(self):
        self.temp.cleanup()

    def test_cli_fetches_official_same_source_record_and_emits_only_version(self):
        self.library.write_bytes(library_bytes("v2.0.0", "v2.0.0-alpha42", source_id=self.source_id))
        output, errors = io.StringIO(), io.StringIO()
        with patch.object(version.urllib.request, "urlopen", return_value=io.BytesIO(b"v2.0.0-alpha42\n")) as fetch:
            with contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors):
                result = version.main(["--library", str(self.library), "--source-dir", str(self.source), "--requested-version", "v2.0.0"])
        self.assertEqual(result, 0, errors.getvalue())
        self.assertEqual(output.getvalue(), "v2.0.0-alpha42\n")
        request = fetch.call_args.args[0]
        self.assertEqual(request.full_url, f"https://duckdb-staging.duckdb.org/{self.source_id}/latest_alpha_version.txt")
        self.assertEqual(fetch.call_args.kwargs["timeout"], 15)

    def test_stable_cli_has_no_network_and_wrong_source_prints_no_flags(self):
        self.library.write_bytes(library_bytes("v2.0.0", source_id=self.source_id))
        with patch.object(version.urllib.request, "urlopen", side_effect=AssertionError("stable must not fetch")):
            self.assertEqual(version.resolve_version(self.library, self.source, "v2.0.0"), "v2.0.0")
        self.library.write_bytes(library_bytes("v2.0.0", source_id="abcdef0123"))
        output, errors = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors):
            result = version.main(["--library", str(self.library), "--source-dir", str(self.source), "--requested-version", "v2.0.0"])
        self.assertEqual(result, 1)
        self.assertEqual(output.getvalue(), "")
        self.assertIn("Source ID", errors.getvalue())

    def test_fetch_failure_is_fail_closed(self):
        self.library.write_bytes(library_bytes("v2.0.0-alpha42", source_id=self.source_id))
        with patch.object(version.urllib.request, "urlopen", side_effect=OSError("offline")):
            with self.assertRaisesRegex(ValueError, "official.*record"):
                version.resolve_version(self.library, self.source, "v2.0.0")

    def make_fixture(self):
        repository = Path(__file__).resolve().parents[1]
        shutil.copyfile(repository / "Makefile", self.root / "Makefile")
        (self.root / "scripts").mkdir()
        shutil.copyfile(repository / "scripts/prebuilt_duckdb_version.py", self.root / "scripts/prebuilt_duckdb_version.py")
        include = self.root / "extension-ci-tools/makefiles"
        include.mkdir(parents=True)
        (include / "duckdb_extension.Makefile").write_text(
            "DUCKDB_SRCDIR ?= ./duckdb\nEXT_RELEASE_FLAGS := -DUPSTREAM_FLAG=1\nrelease:\n\t@python3 -c 'import json,sys; print(json.dumps(sys.argv[1:]))' $(EXT_RELEASE_FLAGS)\n"
        )

    def make(self, prebuilt=True):
        env = os.environ.copy()
        env["DUCKDB_VERSION"] = "v2.0.0"
        env["DUCKDB_SRCDIR"] = str(self.source)
        env["DUCKDB_PREBUILT_LIBRARY"] = str(self.library) if prebuilt else ""
        return subprocess.run(["make", "--no-print-directory", "release"], cwd=self.root, env=env, text=True, capture_output=True)

    def test_make_appends_formal_release_flag_after_include_with_literal_paths(self):
        self.make_fixture()
        self.library.write_bytes(library_bytes("v2.0.0", source_id=self.source_id))
        result = self.make()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout), ["-DUPSTREAM_FLAG=1", "-DDUCKDB_EXPLICIT_VERSION=v2.0.0"])

    def test_make_stops_on_invalid_prebuilt_and_no_prebuilt_is_unchanged(self):
        self.make_fixture()
        self.library.write_bytes(library_bytes("v2.0.0", source_id="abcdef0123"))
        result = self.make()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("prebuilt", result.stderr.lower())
        self.assertEqual(result.stdout, "")
        result = self.make(prebuilt=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout), ["-DUPSTREAM_FLAG=1"])


if __name__ == "__main__":
    unittest.main()
