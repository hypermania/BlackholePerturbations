# Boost 1.84.0 provenance

`external/boost/` is an unchanged copy of the complete public-header directory
`boost_1_84_0/boost/` from the official Boost 1.84.0 source archive. It is not
a `bcp` subset and does not depend on a system Boost installation.

- Release: Boost 1.84.0
- Release date: 2023-12-13
- Archive: `https://archives.boost.io/release/1.84.0/source/boost_1_84_0.tar.bz2`
- Archive metadata: `https://archives.boost.io/release/1.84.0/source/boost_1_84_0.tar.bz2.json`
- Archive SHA-256: `cc4b893acf645c9d4b698e9a0f08ca8846aa5d6c68275c14c3e7949c24109454`
- Upstream commit recorded by the release metadata: `ad09f667e61e18f5c31590941e748ac38e5a81bf`
- Public-header files: 15,689
- Public-header bytes: 144,650,494
- Deterministic public-header tree SHA-256: `84b8450599fa2d3ac0cd8d77962d4827a503c05503990565b500a36c3888f3ee`

The archive also contains an empty `boost/headers/` directory. Git does not
represent empty directories; no file from the official public-header tree is
omitted. The official license is stored alongside the tree as
`external/boost.LICENSE_1_0.txt`, leaving `external/boost/` unchanged.

## Verification

Run:

```bash
make check-boost-vendor
```

The verifier checks the Boost version, file count, total byte count, and a
deterministic digest over every relative path and file content. It then asks
the selected compiler for the complete dependency graph of every Boost entry
point used by the CPU code and tests, failing if any Boost header resolves
outside `external/boost/`.

The tree digest is computed by sorting all relative file paths, writing the
same records as `sha256sum`:

```text
<file SHA-256><two spaces>./<relative path><newline>
```

and taking the SHA-256 of the concatenated records.

## Updating Boost

Update Boost in a separate change from numerical-code changes. Download the
new official release archive and its JSON metadata, verify the archive
SHA-256 before extraction, replace the complete `boost/` directory, update
the constants in `script/verify_boost_vendor.py`, and run the full normal,
production, sanitizer, and CPU application test matrices. Do not regenerate
this dependency with `bcp --scan`: macro-expanded includes are not reliably
discovered, which caused the former incomplete tree.
