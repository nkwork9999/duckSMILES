"""Check translation keys and interpolation placeholders without Xcode."""
import json
import re
from pathlib import Path

resources = Path(__file__).resolve().parents[1] / "Resources"
languages = ["en", "ja", "id", "hi", "de", "fr", "it"]
reference = None
for language in languages:
    source = resources / f"{language}.lproj" / "Localizable.strings"
    pairs = {}
    for line in source.read_text().splitlines():
        if not line.strip():
            continue
        match = re.fullmatch(r'("(?:[^"\\]|\\.)*") = ("(?:[^"\\]|\\.)*");', line)
        assert match, (language, line)
        key, value = (json.loads(match[i]) for i in (1, 2))
        assert key not in pairs, (language, "duplicate", key)
        assert value, (language, "empty", key)
        assert re.findall(r'%(?:lld|@)', key) == re.findall(r'%(?:lld|@)', value), (language, key)
        pairs[key] = value
    if reference is None:
        reference = set(pairs)
    assert set(pairs) == reference, (language, "different keys")
print(json.dumps({"languages": languages, "keys_per_language": len(reference), "matching_keys": True, "matching_placeholders": True}, indent=2))
