"""Stage ``qk``: build kinship (K) and PCA (Q) matrices.

Design notes
------------
* **Deterministic chunk sampling.** The old ``shuf`` made every run different.
  Here the subset of bed chunks is chosen with ``random.Random(seed)``, so runs
  are reproducible.
* **Order safety.** All bed chunks share the canonical fam order, so the merged
  ``QK_bed`` is produced by *concatenating* chunk beds (see :func:`bed.concat_beds`)
  rather than ``plink --merge-list`` -- which would silently reorder samples.
  The PCA ``.eigenvec`` IDs are then checked against the canonical order.
* External K/PCA is recommended for real analyses; if supplied we verify the PCA
  order when it carries IDs, and otherwise document the requirement.
"""

from __future__ import annotations

import glob
import math
import os
import random
import re
import shutil

from ..bed import concat_beds, fam_sample_ids
from ..context import RunContext
from ..errors import (
    AkgwasError,
    E_GEMMA_FAILED,
    E_INPUT_MISSING,
    E_PLINK_FAILED,
    E_SAMPLE_ORDER_MISMATCH,
    E_SAMPLE_ORDER_UNVERIFIABLE,
)
from ..tools import atomic_write_text, ensure_dir, remove_if_exists, require_tools

STAGE = "qk"

KINSHIP_NAME = "QK_kinship.txt"
PCA_NAME = "QK_pca1-3.txt"
SAMPLES_NAME = "QK_samples.txt"


def run(ctx: RunContext) -> None:
    cfg = ctx.cfg
    qk_dir = cfg.qk_path
    bed_dir = cfg.bed_path
    inputs = [bed_dir] if not cfg.qk.use_external else [cfg.external_kinship_path, cfg.external_pca_path]

    if not ctx.force and ctx.markers.is_done(STAGE, inputs):
        ctx.log("qk: marker valid -> skip")
        return
    ctx.markers.clear(STAGE)
    ensure_dir(qk_dir)

    canonical = _canonical_samples(bed_dir)
    if not canonical:
        raise AkgwasError(E_INPUT_MISSING, f"no bed chunks found under {bed_dir}")
    atomic_write_text(os.path.join(qk_dir, SAMPLES_NAME), "\n".join(canonical) + "\n")

    if cfg.qk.use_external:
        _setup_external(ctx, canonical)
    else:
        _build_internal(ctx, canonical)

    if not ctx.dry_run:
        ctx.markers.mark_done(STAGE, inputs)


def _canonical_samples(bed_dir: str) -> list[str]:
    order_file = os.path.join(bed_dir, "sample_id_inorder.txt")
    if os.path.exists(order_file):
        with open(order_file, encoding="utf-8") as fh:
            return [ln.strip() for ln in fh if ln.strip()]
    chunks = sorted(glob.glob(os.path.join(bed_dir, "matrix_*.bed")))
    if chunks:
        return fam_sample_ids(chunks[0][:-4])
    return []


def _chunks(bed_dir: str) -> list[str]:
    return sorted(p[:-4] for p in glob.glob(os.path.join(bed_dir, "matrix_*.bed")))


# ---------------------------------------------------------------------------
# Internal (k-mer) K/PCA
# ---------------------------------------------------------------------------

def _select_chunks(chunks: list[str], fraction: float, seed: int) -> list[str]:
    n = len(chunks)
    k = int(math.ceil(n * fraction))
    k = max(1, min(k, n))
    rng = random.Random(seed)
    picked = sorted(rng.sample(range(n), k))
    return [chunks[i] for i in picked]


def _build_internal(ctx: RunContext, canonical: list[str]) -> None:
    cfg = ctx.cfg
    qk_dir = cfg.qk_path
    bed_dir = cfg.bed_path
    require_tools(["plink", "gemma"])

    chunks = _chunks(bed_dir)
    if not chunks:
        raise AkgwasError(E_INPUT_MISSING, f"no bed chunks found under {bed_dir}")
    selected = _select_chunks(chunks, cfg.qk.sample_fraction, cfg.qk.seed)
    ctx.log(f"qk: using {len(selected)}/{len(chunks)} chunks (seed={cfg.qk.seed})")

    qk_bed = os.path.join(qk_dir, "QK_bed")
    concat_beds(selected, qk_bed)
    _patch_fam_phenotype(qk_bed + ".fam", value="1")

    # Verify concatenation preserved the canonical order.
    merged_ids = fam_sample_ids(qk_bed + ".fam")
    if merged_ids != canonical:
        raise AkgwasError(
            E_SAMPLE_ORDER_MISMATCH,
            "internal QK_bed sample order does not match the canonical bed order",
        )

    # PCA via plink.
    ctx.run_tool(
        ["plink", "--bfile", qk_bed, "--pca", "3", "--threads", cfg.qk.threads,
         "--out", os.path.join(qk_dir, "QK_pca")],
        error_code=E_PLINK_FAILED,
    )
    _verify_and_write_pca(os.path.join(qk_dir, "QK_pca.eigenvec"), os.path.join(qk_dir, PCA_NAME), canonical)

    # Kinship via gemma -gk.
    ctx.run_tool(
        ["gemma", "-bfile", qk_bed, "-gk", cfg.qk.gk, "-miss", "1.0", "-maf", "0.0",
         "-r2", "1.0", "-o", "QK_kinship"],
        cwd=qk_dir,
        error_code=E_GEMMA_FAILED,
    )
    _collect_gemma_outputs(qk_dir, "QK_kinship", KINSHIP_NAME, cfg.qk.gk)

    # QK_bed is an intermediate; drop the heavy files, keep the arrays.
    for ext in ("bed", "bim", "fam", "nosex", "log"):
        remove_if_exists(f"{qk_bed}.{ext}")


def _patch_fam_phenotype(fam_path: str, value: str) -> None:
    """Set the 6th (phenotype) column; GEMMA is happier with a non-missing value."""
    with open(fam_path, encoding="utf-8") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    out = []
    for line in lines:
        parts = line.split()
        while len(parts) < 6:
            parts.append("0")
        parts[5] = value
        out.append(" ".join(parts))
    atomic_write_text(fam_path, "\n".join(out) + "\n")


def _collect_gemma_outputs(qk_dir: str, prefix: str, target_name: str, gk: int) -> None:
    """Move GEMMA's relatedness matrix up and rename it to a stable name.

    GEMMA names the matrix after the estimator: ``<prefix>.cXX.txt`` for
    ``-gk 1`` (centred, the default) and ``<prefix>.sXX.txt`` for ``-gk 2``
    (standardized). ``XX`` is a *literal* part of the name -- GEMMA does **not**
    substitute the sample count, so the real file is ``QK_kinship.cXX.txt``, not
    ``QK_kinship.c20.txt``. Only those two shapes are accepted -- in particular
    the run log (``<prefix>.log.txt``) must never be mistaken for a matrix. The
    old code globbed ``<prefix>.*.txt`` and took the lexically-first hit across
    both ``output/`` and ``qk_dir``; on a second run the leftover
    ``qk_dir/<prefix>.log.txt`` ('Q' < 'o') won and got copied over
    ``QK_kinship.txt``, so GEMMA later rejected it ("number of columns in the
    kinship file does not match the number of individuals").
    """
    kind = "s" if gk == 2 else "c"
    # ``XX`` is literal in the filename; ``\d+`` is tolerated for safety across
    # GEMMA builds, but the real GEMMA output is ``<prefix>.cXX.txt``.
    matrix_re = re.compile(rf"^{re.escape(prefix)}\.{kind}(?:XX|\d+)\.txt$")
    candidates = [
        path
        for folder in (os.path.join(qk_dir, "output"), qk_dir)
        for path in glob.glob(os.path.join(folder, f"{prefix}.*.txt"))
        if matrix_re.match(os.path.basename(path))
    ]
    if not candidates:
        raise AkgwasError(
            E_GEMMA_FAILED,
            f"gemma -gk {gk} did not produce a {prefix}.{kind}XX.txt kinship "
            f"matrix under {qk_dir}",
        )
    shutil.move(sorted(candidates)[0], os.path.join(qk_dir, target_name))
    # Tidy the gemma output/ dir.
    outdir = os.path.join(qk_dir, "output")
    if os.path.isdir(outdir):
        for name in os.listdir(outdir):
            shutil.move(os.path.join(outdir, name), os.path.join(qk_dir, name))
        os.rmdir(outdir)


def _verify_and_write_pca(eigenvec: str, dest: str, canonical: list[str]) -> None:
    if not os.path.exists(eigenvec):
        raise AkgwasError(E_PLINK_FAILED, f"plink --pca did not write {eigenvec}")
    rows = []
    with open(eigenvec, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            rows.append(parts)
    ids = [r[1] if len(r) > 1 else r[0] for r in rows]
    if ids != canonical:
        raise AkgwasError(
            E_SAMPLE_ORDER_MISMATCH,
            "PCA eigenvector sample order does not match the bed fam order",
        )
    body = "\n".join(" ".join(r[2:5]) for r in rows if len(r) >= 5) + "\n"
    atomic_write_text(dest, body)


# ---------------------------------------------------------------------------
# External K/PCA
# ---------------------------------------------------------------------------

def _setup_external(ctx: RunContext, canonical: list[str]) -> None:
    cfg = ctx.cfg
    qk_dir = cfg.qk_path
    kinship_src = cfg.external_kinship_path
    pca_src = cfg.external_pca_path
    if not os.path.exists(kinship_src):
        raise AkgwasError(E_INPUT_MISSING, f"external kinship not found: {kinship_src}")
    if not os.path.exists(pca_src):
        raise AkgwasError(E_INPUT_MISSING, f"external PCA not found: {pca_src}")

    shutil.copyfile(kinship_src, os.path.join(qk_dir, KINSHIP_NAME))

    if cfg.qk.pca_has_id:
        _import_external_pca_with_ids(pca_src, os.path.join(qk_dir, PCA_NAME), canonical)
    else:
        shutil.copyfile(pca_src, os.path.join(qk_dir, PCA_NAME))
        n_rows = sum(1 for ln in open(pca_src, encoding="utf-8") if ln.strip())
        if n_rows != len(canonical):
            raise AkgwasError(
                E_SAMPLE_ORDER_UNVERIFIABLE,
                f"external PCA has {n_rows} rows but the bed has {len(canonical)} samples; "
                "provide PCA with sample IDs and pass --pca-has-id, or fix the file",
            )
        ctx.logger.warn(
            "external PCA has no sample IDs; assuming its row order already matches the bed fam order"
        )

    ctx.log("qk: using external K/PCA; kinship row order is assumed as documented in README")


def _import_external_pca_with_ids(src: str, dest: str, canonical: list[str]) -> None:
    rows = {}
    header_skipped = False
    with open(src, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.split()
            if not header_skipped and not _is_numeric(parts[2] if len(parts) > 2 else ""):
                header_skipped = True
                continue
            if len(parts) < 3:
                continue
            rows[parts[1]] = parts[2:5]
    missing = [s for s in canonical if s not in rows]
    if missing:
        raise AkgwasError(
            E_SAMPLE_ORDER_UNVERIFIABLE,
            f"external PCA is missing {len(missing)} bed samples, e.g. {missing[:3]}",
        )
    body = "\n".join(" ".join(rows[s]) for s in canonical) + "\n"
    atomic_write_text(dest, body)


def _is_numeric(text: str) -> bool:
    try:
        float(text)
        return True
    except ValueError:
        return False
