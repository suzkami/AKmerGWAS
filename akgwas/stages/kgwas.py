"""Stage ``kgwas``: per (chunk x trait) GEMMA association + significant-k-mer extraction.

Before GWAS runs, the multi-column phenotype file is validated against the bed
``.fam`` ID set, reordered to the fam order, and split into one single-column
file per trait (see :mod:`akgwas.pheno`). Each (chunk, trait) pair is an
independent single-threaded GEMMA call; ``--cores`` fans them out.
"""

from __future__ import annotations

import glob
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

from ..bed import fam_sample_ids
from ..context import RunContext
from ..errors import AkgwasError, E_ASSOC_MALFORMED, E_GEMMA_FAILED, E_INPUT_MISSING
from ..pheno import export_per_trait
from ..tools import atomic_write_text, ensure_dir, require_tools

STAGE = "kgwas"

# .assoc.txt column indices (0-based): rs, af, beta, se, p_wald
_ASSOC_COLS = [1, 6, 7, 8, 11]


def run(ctx: RunContext) -> None:
    cfg = ctx.cfg
    bed_dir = cfg.bed_path
    qk_dir = cfg.qk_path
    temp_dir = cfg.temp_path
    pheno_single_dir = cfg.pheno_single_path

    kinship = os.path.join(qk_dir, "QK_kinship.txt")
    pca = os.path.join(qk_dir, "QK_pca1-3.txt")
    inputs = [bed_dir, cfg.pheno_path, kinship, pca]

    key = ctx.marker_key(STAGE)
    if not ctx.force and ctx.markers.is_done(key, inputs):
        ctx.log("kgwas: marker valid -> skip")
        return
    ctx.markers.clear(key)

    for path in (kinship, pca):
        if not os.path.exists(path):
            raise AkgwasError(E_INPUT_MISSING, f"missing K/PCA file: {path} (run the qk stage first)")

    chunks = sorted(p[:-4] for p in glob.glob(os.path.join(bed_dir, "matrix_*.bed")))
    if not chunks:
        raise AkgwasError(E_INPUT_MISSING, f"no bed chunks under {bed_dir}")

    # Validate + reorder + export the phenotype columns.
    fam_path = chunks[0] + ".fam"
    traits, trait_files = export_per_trait(cfg.pheno_path, fam_path, pheno_single_dir)
    ctx.log(f"kgwas: {len(traits)} trait(s) x {len(chunks)} chunk(s)")

    if not ctx.dry_run:
        require_tools(["gemma"])
        for chunk in chunks:
            ensure_dir(os.path.join(temp_dir, os.path.basename(chunk)))

    jobs = [(chunk, trait, trait_files[trait]) for chunk in chunks for trait in traits]

    def work(job):
        chunk, trait, pheno_file = job
        _run_one(ctx, chunk, trait, pheno_file, kinship, pca)

    if ctx.dry_run:
        for chunk, trait, _ in jobs:
            ctx.run_tool(["gemma", "-lmm", "1", "...", "-o", trait])  # logs the shape
        return

    errors = []
    with ThreadPoolExecutor(max_workers=max(1, ctx.cores)) as pool:
        futures = {pool.submit(work, job): job for job in jobs}
        for future in as_completed(futures):
            chunk, trait, _ = futures[future]
            try:
                future.result()
            except Exception as exc:  # noqa: BLE001 - re-raised below
                errors.append(f"{os.path.basename(chunk)}/{trait}: {exc}")
    if errors:
        raise AkgwasError(
            E_ASSOC_MALFORMED, "kgwas failures:\n  " + "\n  ".join(errors[:10])
        )

    ctx.markers.mark_done(key, inputs)
    ctx.log("kgwas: done")


def _run_one(ctx: RunContext, chunk: str, trait: str, pheno_file: str, kinship: str, pca: str) -> None:
    cfg = ctx.cfg
    work_dir = os.path.join(cfg.temp_path, os.path.basename(chunk))
    ensure_dir(work_dir)

    env = {
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    try:
        ctx.run_tool(
            ["gemma", "-bfile", chunk, "-k", kinship, "-p", pheno_file, "-c", pca,
             "-lmm", "1", "-o", trait],
            cwd=work_dir,
            env=env,
            error_code=E_GEMMA_FAILED,
        )
    except AkgwasError as exc:
        # GEMMA can exit non-zero for numerical reasons (e.g. a rank-deficient
        # GRM on tiny datasets). Surface the tail of its log so the user sees the
        # actual GSL/Brent message rather than a bare "exited N".
        raise AkgwasError(E_GEMMA_FAILED, exc.message + _gemma_log_tail(work_dir, trait)) from exc

    assoc = os.path.join(work_dir, "output", f"{trait}.assoc.txt")
    if not os.path.exists(assoc):
        raise AkgwasError(E_ASSOC_MALFORMED, f"gemma did not produce {assoc}")

    sig_rows = _significant_rows(assoc, cfg.kgwas.threshold)

    sigk_dir = os.path.join(work_dir, "sigk")
    fasta_dir = os.path.join(work_dir, "fasta")
    ensure_dir(sigk_dir)
    ensure_dir(fasta_dir)

    # sigK columns: k-mer, af, beta, se, p_wald (no header).
    atomic_write_text(
        os.path.join(sigk_dir, f"{trait}_sigK.txt"),
        "".join(row + "\n" for row in sig_rows),
    )
    kmers = [row.split("\t")[0] for row in sig_rows]
    atomic_write_text(
        os.path.join(fasta_dir, f"{trait}_fasta.txt"),
        "".join(f">{k}\n{k}\n" for k in kmers),
    )


def _gemma_log_tail(work_dir: str, trait: str, n: int = 10) -> str:
    """Return the last ``n`` non-empty lines of GEMMA's log, or "" if absent."""
    log = os.path.join(work_dir, "output", f"{trait}.log.txt")
    if not os.path.exists(log):
        return ""
    try:
        with open(log, encoding="utf-8") as fh:
            lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    except OSError:
        return ""
    if not lines:
        return ""
    return "\n    gemma log tail:\n      " + "\n      ".join(lines[-n:])


def _significant_rows(assoc_path: str, threshold: float) -> list[str]:
    rows = []
    with open(assoc_path, encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            if i == 0 or not line.strip():
                continue
            if len(parts) <= max(_ASSOC_COLS):
                raise AkgwasError(
                    E_ASSOC_MALFORMED,
                    f"{assoc_path}: expected >= {max(_ASSOC_COLS) + 1} columns, found {len(parts)}",
                )
            try:
                p = float(parts[_ASSOC_COLS[-1]])
            except ValueError:
                continue
            if p <= threshold:
                rows.append("\t".join(parts[c] for c in _ASSOC_COLS))
    rows.sort(key=lambda r: r.split("\t")[0])
    return rows
