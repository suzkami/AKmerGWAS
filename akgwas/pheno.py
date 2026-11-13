"""Phenotype handling for the v1.0.0 multi-column format.

The phenotype file is now a **single table**:

    ID        ft10     ft16
    A88       61.25    22.1
    A108      70       NA
    ...

* the first column holds sample IDs,
* the header row holds phenotype (trait) names,
* missing values are written as ``NA``.

GWAS requires one single-column, header-less file per trait. We therefore
(1) verify the sample-ID *set* equals the bed ``.fam`` set, (2) reorder rows to
the ``.fam`` order (the single canonical order for the whole run), and
(3) export each trait column to ``<run_root>/pheno_single/<trait>.txt``.

This is the *only* place sample rows are reordered, so "row alignment" reduces
to "conform everything to the bed fam order".
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass

from .bed import fam_sample_ids
from .errors import AkgwasError, E_PHENO_MALFORMED, E_SAMPLE_ID_MISMATCH
from .tools import atomic_write_text, ensure_dir

_SAFE = re.compile(r"[^A-Za-z0-9._-]+")


def sanitize_trait(name: str) -> str:
    """Make a trait name safe to use as a file name stem."""
    cleaned = _SAFE.sub("_", name.strip())
    return cleaned or "trait"


@dataclass
class PhenoTable:
    id_header: str
    traits: list[str]
    ids: list[str]
    values: dict[str, list[str]]  # trait -> values aligned with ``ids``


def read_pheno(path: str) -> PhenoTable:
    # Keep original line numbers for error messages; blank lines and ``#``
    # comments are ignored (consistent with the fof / reference-list readers).
    rows: list[tuple[int, str]] = []
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, start=1):
            text = line.rstrip("\n")
            if text.strip() and not text.lstrip().startswith("#"):
                rows.append((lineno, text))
    if not rows:
        raise AkgwasError(E_PHENO_MALFORMED, f"phenotype file is empty: {path}")

    header = rows[0][1].split()
    if len(header) < 2:
        raise AkgwasError(
            E_PHENO_MALFORMED,
            "phenotype header must have an ID column plus at least one trait column "
            "(first column = sample ID, remaining columns = trait names)",
        )
    id_header, traits = header[0], [sanitize_trait(t) for t in header[1:]]

    ids: list[str] = []
    columns: dict[str, list[str]] = {t: [] for t in traits}
    seen: set[str] = set()
    for lineno, line in rows[1:]:
        parts = line.split()
        if len(parts) != len(header):
            raise AkgwasError(
                E_PHENO_MALFORMED,
                f"{path}:{lineno}: expected {len(header)} columns, found {len(parts)}",
            )
        sample = parts[0]
        if sample in seen:
            raise AkgwasError(E_PHENO_MALFORMED, f"{path}:{lineno}: duplicate sample id {sample!r}")
        seen.add(sample)
        ids.append(sample)
        for trait, value in zip(traits, parts[1:]):
            columns[trait].append(value)

    return PhenoTable(id_header=id_header, traits=traits, ids=ids, values=columns)


def export_per_trait(pheno_path: str, fam_path: str, outdir: str) -> tuple[list[str], dict[str, str]]:
    """Validate + reorder + export one single-column file per trait.

    Returns ``(trait_names, {trait: path})``. Raises on an ID-set mismatch.
    """
    table = read_pheno(pheno_path)
    fam_ids = fam_sample_ids(fam_path)

    fam_set = set(fam_ids)
    pheno_set = set(table.ids)
    if fam_set != pheno_set:
        missing = sorted(fam_set - pheno_set)
        extra = sorted(pheno_set - fam_set)
        detail = []
        if missing:
            detail.append(f"missing from phenotype ({len(missing)}): {', '.join(missing[:5])}...")
        if extra:
            detail.append(f"extra in phenotype ({len(extra)}): {', '.join(extra[:5])}...")
        raise AkgwasError(
            E_SAMPLE_ID_MISMATCH,
            "sample IDs in the phenotype file differ from the bed .fam. " + "; ".join(detail),
        )

    row_index = {sid: i for i, sid in enumerate(table.ids)}
    order = [row_index[sid] for sid in fam_ids]

    ensure_dir(outdir)
    paths: dict[str, str] = {}
    for trait in table.traits:
        col = table.values[trait]
        body = "\n".join(col[i] for i in order) + "\n"
        dest = os.path.join(outdir, f"{trait}.txt")
        atomic_write_text(dest, body)
        paths[trait] = dest
    return table.traits, paths
