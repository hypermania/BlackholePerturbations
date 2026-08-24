#!/usr/bin/env python3
"""Verify the vendored Boost 1.84 public-header tree and include closure."""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
from pathlib import Path


BOOST_VERSION = 108400
EXPECTED_FILE_COUNT = 15_689
EXPECTED_TOTAL_BYTES = 144_650_494
EXPECTED_TREE_SHA256 = (
    "84b8450599fa2d3ac0cd8d77962d4827a503c05503990565b500a36c3888f3ee"
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
BOOST_ROOT = PROJECT_ROOT / "external" / "boost"

# Every Boost entry point included directly by the current CPU code and tests.
PROJECT_BOOST_HEADERS = (
    "boost/lockfree/queue.hpp",
    "boost/math/constants/constants.hpp",
    "boost/math/special_functions/lambert_w.hpp",
    "boost/math/tools/roots.hpp",
    "boost/multiprecision/cpp_bin_float.hpp",
    "boost/multiprecision/eigen.hpp",
    "boost/multiprecision/float128.hpp",
    "boost/numeric/odeint.hpp",
    "boost/numeric/odeint/algebra/algebra_dispatcher.hpp",
    "boost/numeric/odeint/algebra/operations_dispatcher.hpp",
    "boost/numeric/odeint/external/eigen/eigen.hpp",
    "boost/numeric/odeint/util/copy.hpp",
    "boost/numeric/odeint/util/resize.hpp",
    "boost/numeric/odeint/util/same_size.hpp",
    "boost/pfr.hpp",
    "boost/range.hpp",
    "boost/type_index.hpp",
    # Regression check for the header missed by bcp --scan.
    "boost/typeof/incr_registration_group.hpp",
)


def tree_fingerprint(root: Path) -> tuple[int, int, str]:
    """Return file count, byte count, and a path-and-content tree digest."""
    relative_paths = sorted(
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file()
    )
    digest = hashlib.sha256()
    total_bytes = 0

    for relative_path in relative_paths:
        path = root / relative_path
        contents = path.read_bytes()
        total_bytes += len(contents)
        file_digest = hashlib.sha256(contents).hexdigest()
        # This matches the deterministic shell construction documented in the
        # provenance file: sorted `sha256sum` records hashed once more.
        digest.update(f"{file_digest}  ./{relative_path}\n".encode())

    return len(relative_paths), total_bytes, digest.hexdigest()


def verify_tree() -> None:
    if not BOOST_ROOT.is_dir():
        raise RuntimeError(f"vendored Boost directory is missing: {BOOST_ROOT}")

    non_files = [
        path
        for path in BOOST_ROOT.rglob("*")
        if path.is_symlink() or (not path.is_file() and not path.is_dir())
    ]
    if non_files:
        raise RuntimeError(f"unexpected non-file entries under Boost: {non_files[:5]}")

    version_text = (BOOST_ROOT / "version.hpp").read_text()
    version_declaration = f"#define BOOST_VERSION {BOOST_VERSION}"
    if version_declaration not in version_text:
        raise RuntimeError(f"expected {version_declaration!r} in boost/version.hpp")

    file_count, total_bytes, tree_sha256 = tree_fingerprint(BOOST_ROOT)
    actual = (file_count, total_bytes, tree_sha256)
    expected = (
        EXPECTED_FILE_COUNT,
        EXPECTED_TOTAL_BYTES,
        EXPECTED_TREE_SHA256,
    )
    if actual != expected:
        raise RuntimeError(
            "vendored Boost tree differs from the official 1.84.0 public "
            "headers:\n"
            f"  expected files/bytes/SHA-256: {expected}\n"
            f"  actual files/bytes/SHA-256:   {actual}"
        )

    print(
        "Boost 1.84.0 tree verified: "
        f"{file_count} files, {total_bytes} bytes, {tree_sha256}"
    )


def is_within(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def verify_include_closure(compiler: str) -> None:
    missing_headers = [
        header
        for header in PROJECT_BOOST_HEADERS
        if not (PROJECT_ROOT / "external" / header).is_file()
    ]
    if missing_headers:
        raise RuntimeError(f"project Boost entry points are missing: {missing_headers}")

    source = "".join(f"#include <{header}>\n" for header in PROJECT_BOOST_HEADERS)
    command = [
        compiler,
        "-std=c++20",
        "-Iexternal",
        "-M",
        "-x",
        "c++",
        "-",
    ]
    result = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        input=source,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "Boost dependency preprocessing failed:\n"
            f"command: {' '.join(command)}\n{result.stderr}"
        )

    dependency_text = result.stdout.replace("\\\n", " ")
    boost_dependencies: set[Path] = set()
    for token in dependency_text.split():
        token = token.rstrip("\\")
        candidate = Path(token)
        if "boost" not in candidate.parts:
            continue
        resolved = candidate.resolve() if candidate.is_absolute() else (PROJECT_ROOT / candidate).resolve()
        boost_dependencies.add(resolved)

    boost_root = BOOST_ROOT.resolve()
    leaked_headers = sorted(
        path for path in boost_dependencies if not is_within(path, boost_root)
    )
    if leaked_headers:
        formatted = "\n".join(f"  {path}" for path in leaked_headers)
        raise RuntimeError(f"Boost headers resolved outside external/boost:\n{formatted}")

    print(
        "Boost include closure verified: "
        f"{len(boost_dependencies)} transitive headers all resolve under external/boost"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compiler",
        default="g++",
        help="C++ compiler used for the transitive include audit (default: g++)",
    )
    args = parser.parse_args()

    try:
        verify_tree()
        verify_include_closure(args.compiler)
    except RuntimeError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
