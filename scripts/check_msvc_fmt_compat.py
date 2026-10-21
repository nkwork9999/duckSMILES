#!/usr/bin/env python3
"""Compile a C translation unit with the MSVC fmt shim force-included.

DuckDB 2 generates C entry points in addition to C++ extension sources. The
shim must not include C++ STL headers or change C preprocessor state there.
"""

import os
from pathlib import Path
import shlex
import subprocess
import tempfile


def main():
    header = Path(__file__).resolve().parents[1] / "cmake/msvc_fmt_compat.hpp"
    compiler = shlex.split(os.environ.get("CC", "cc"))
    with tempfile.TemporaryDirectory(prefix="ducksmiles-msvc-c-") as temp:
        source = Path(temp) / "generated_entry.c"
        source.write_text(
            "#if !defined(_SECURE_SCL) || _SECURE_SCL != 1\n"
            '#error "The C translation unit must retain its macro state"\n'
            "#endif\n"
            "int generated_entry(void) { return 0; }\n"
        )
        subprocess.run(
            compiler
            + ["-D_MSC_VER=1950", "-D_SECURE_SCL=1", "-include", str(header),
               "-fsyntax-only", str(source)],
            check=True,
        )
    print("PASS: MSVC fmt shim is safe in generated C translation units")


if __name__ == "__main__":
    main()
