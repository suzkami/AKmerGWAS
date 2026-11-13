"""Parsers for the plain-text input lists (fof, reference list)."""

from __future__ import annotations

from dataclasses import dataclass

from .errors import AkgwasError, E_FOF_MALFORMED, E_REF_LIST_EMPTY


@dataclass
class FofEntry:
    name: str
    files: list[str]


def read_fof(path: str) -> list[FofEntry]:
    """Parse ``fof.txt``: ``name : file1 ; file2 ; ...`` (one sample per line)."""
    entries: list[FofEntry] = []
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" not in line:
                raise AkgwasError(
                    E_FOF_MALFORMED,
                    f"{path}:{lineno}: expected `name : file1 ; file2`, got {line!r}",
                )
            name, rest = line.split(":", 1)
            name = name.strip()
            files = [x.strip() for x in rest.split(";") if x.strip()]
            if not name or not files:
                raise AkgwasError(
                    E_FOF_MALFORMED, f"{path}:{lineno}: empty sample name or file list"
                )
            entries.append(FofEntry(name=name, files=files))
    if not entries:
        raise AkgwasError(E_FOF_MALFORMED, f"{path}: no samples found")
    return entries


def read_fof_samples(path: str) -> list[str]:
    return [e.name for e in read_fof(path)]


def read_ref_list(path: str) -> list[str]:
    """Read reference-index entries (one per line, ``#`` comments ignored)."""
    refs: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                refs.append(line)
    if not refs:
        raise AkgwasError(E_REF_LIST_EMPTY, f"{path}: no reference indices listed")
    return refs
