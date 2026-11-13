"""Stage ``collect``: merge per-chunk GEMMA results into per-trait outputs.

Fixes the original collector's P0 bug (append-never-truncate + looping over all
traits while the scheduler also parallelised over traits). Here collection is a
single, sequential, per-trait pass that **writes fresh** (atomic replace).

The per-trait ``siggeno`` bed is built by pulling each chunk's significant
variants out of that chunk bed and concatenating -- no ``plink --merge-list``,
so the fam order is preserved by construction.

The full, all-k-mers association dump (``gemma_combine/<trait>.assoc_combined.txt``)
is the biggest output of a run, so it is **opt-in**: it is only written when
``--gemma-combine`` is passed (``cfg.collect.gemma_combine``). When off we also
skip *reading* every chunk's ``.assoc.txt``, which is itself a large I/O cost.
"""

from __future__ import annotations

import glob
import os

from ..bed import concat_beds, subset_bed
from ..context import RunContext
from ..errors import AkgwasError, E_COLLECT_EMPTY
from ..tools import atomic_write_text, ensure_dir, remove_if_exists

STAGE = "collect"

_ASSOC_COLS = [1, 6, 7, 8, 11]


def run(ctx: RunContext) -> None:
    cfg = ctx.cfg
    temp_dir = cfg.temp_path
    final_dir = cfg.final_path
    pheno_single = cfg.pheno_single_path
    inputs = [temp_dir, pheno_single]

    key = ctx.marker_key(STAGE)
    if not ctx.force and ctx.markers.is_done(key, inputs):
        ctx.log("collect: marker valid -> skip")
        return
    ctx.markers.clear(key)

    chunks = sorted(
        os.path.basename(p)
        for p in glob.glob(os.path.join(temp_dir, "*"))
        if os.path.isdir(p)
    )
    traits = sorted(
        os.path.splitext(os.path.basename(p))[0]
        for p in glob.glob(os.path.join(pheno_single, "*.txt"))
    )
    if not traits:
        raise AkgwasError(E_COLLECT_EMPTY, f"no per-trait phenotype files under {pheno_single}")

    subs = ["sigk", "fasta", "siggeno", "logs"]
    if cfg.collect.gemma_combine:
        subs.append("gemma_combine")
    for sub in subs:
        ensure_dir(os.path.join(final_dir, sub))

    if ctx.dry_run:
        ctx.log(f"DRY-RUN: collect {len(traits)} trait(s) from {len(chunks)} chunk(s)")
        return

    for trait in traits:
        _collect_trait(ctx, trait, chunks)

    ctx.markers.mark_done(key, inputs)
    ctx.log("collect: done")


def _collect_trait(ctx: RunContext, trait: str, chunks: list[str]) -> None:
    cfg = ctx.cfg
    temp_dir = cfg.temp_path
    final_dir = cfg.final_path
    keep_full = cfg.collect.gemma_combine

    sigk_rows: list[str] = []
    assoc_rows: list[str] = []
    fasta_parts: list[str] = []
    log_parts: list[str] = []

    # Per-chunk significant k-mers (so we can subset the matching chunk bed).
    per_chunk_kmers: dict[str, list[str]] = {}
    for chunk in chunks:
        cdir = os.path.join(temp_dir, chunk)
        sigk_file = os.path.join(cdir, "sigk", f"{trait}_sigK.txt")
        kmers: list[str] = []
        if os.path.exists(sigk_file):
            with open(sigk_file, encoding="utf-8") as fh:
                for line in fh:
                    if line.strip():
                        sigk_rows.append(line.rstrip("\n"))
                        kmers.append(line.split("\t")[0])
        per_chunk_kmers[chunk] = kmers

        assoc_file = os.path.join(cdir, "output", f"{trait}.assoc.txt")
        if keep_full and os.path.exists(assoc_file):
            assoc_rows.extend(_read_assoc(assoc_file))

        fasta_file = os.path.join(cdir, "fasta", f"{trait}_fasta.txt")
        if os.path.exists(fasta_file):
            with open(fasta_file, encoding="utf-8") as fh:
                fasta_parts.append(fh.read())

        log_file = os.path.join(cdir, "output", f"{trait}.log.txt")
        if os.path.exists(log_file):
            with open(log_file, encoding="utf-8") as fh:
                log_parts.append(f"# {cdir}\n" + fh.read())

    sigk_rows.sort(key=lambda r: r.split("\t")[0])
    atomic_write_text(os.path.join(final_dir, "sigk", f"{trait}_sigK.txt"),
                      "".join(r + "\n" for r in sigk_rows))
    if keep_full:
        atomic_write_text(os.path.join(final_dir, "gemma_combine", f"{trait}.assoc_combined.txt"),
                          "".join(r + "\n" for r in assoc_rows))
    atomic_write_text(os.path.join(final_dir, "fasta", f"{trait}_fasta.txt"),
                      "".join(fasta_parts))
    atomic_write_text(os.path.join(final_dir, "logs", f"{trait}.log_combined.txt"),
                      "".join(log_parts))

    _build_siggeno(ctx, trait, per_chunk_kmers)
    ctx.log(f"collect[{trait}]: {len(sigk_rows)} significant k-mer(s)")


def _read_assoc(path: str) -> list[str]:
    rows = []
    with open(path, encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            if i == 0 or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(_ASSOC_COLS):
                continue
            rows.append("\t".join(parts[c] for c in _ASSOC_COLS))
    return rows


def _build_siggeno(ctx: RunContext, trait: str, per_chunk_kmers: dict[str, list[str]]) -> None:
    cfg = ctx.cfg
    bed_dir = cfg.bed_path
    siggeno_dir = os.path.join(cfg.final_path, "siggeno")
    out_prefix = os.path.join(siggeno_dir, f"{trait}_bed")

    subset_prefixes = []
    for chunk, kmers in per_chunk_kmers.items():
        if not kmers:
            continue
        src = os.path.join(bed_dir, chunk)
        if not os.path.exists(src + ".bed"):
            continue
        tmp = os.path.join(siggeno_dir, f".tmp_{trait}_{chunk}")
        if subset_bed(src, kmers, tmp) > 0:
            subset_prefixes.append(tmp)
        else:  # none of this chunk's k-mers were significant -> drop the empty bed
            for ext in ("bed", "bim", "fam"):
                remove_if_exists(f"{tmp}.{ext}")

    if not subset_prefixes:
        ctx.log(f"collect[{trait}]: no significant k-mers -> no siggeno bed")
        return

    subset_prefixes.sort()
    concat_beds(subset_prefixes, out_prefix)
    for prefix in subset_prefixes:
        for ext in ("bed", "bim", "fam"):
            remove_if_exists(f"{prefix}.{ext}")
