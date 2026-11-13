"""Stage ``alignment``: map significant k-mers back to reference genomes.

Loop order (requirement 7): the **outer** unit is the reference, the inner unit
is the trait. Running every trait against reference 1, then reference 2, ...,
keeps each reference's index hot in the page cache across traits; interleaving
references would thrash the cache. References run in parallel (``--cores``);
within one reference the traits run sequentially.

Filtering (requirement 2): keep every full-length exact hit -- including
secondary/supplementary (flag 256/2048) and reverse-strand (flag 16) records --
drop only true unmapped records (``flag & 4``), and drop an entire k-mer if it
has more than ``alignment.max_hits`` positions.
"""

from __future__ import annotations

import glob
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass

from ..context import RunContext
from ..errors import AkgwasError, E_ALIGNER_UNKNOWN, E_ALIGN_FAILED, E_INPUT_MISSING
from ..inputs import read_ref_list
from ..tools import atomic_write_text, ensure_dir, require_tools

STAGE = "alignment"

_SEQ_EXTS = (".fasta.gz", ".fa.gz", ".fasta", ".fa", ".fna.gz", ".fna", ".gz", ".sa")


@dataclass
class RefInfo:
    entry: str
    index: str
    name: str


def run(ctx: RunContext) -> None:
    cfg = ctx.cfg
    final_dir = cfg.final_path
    traits = _traits(final_dir)
    refs = _refs(ctx)
    inputs = [os.path.join(final_dir, "sigk"), os.path.join(final_dir, "fasta"), cfg.ref_list_path]

    key = ctx.marker_key(STAGE)
    if not ctx.force and ctx.markers.is_done(key, inputs):
        ctx.log("alignment: marker valid -> skip")
        return
    ctx.markers.clear(key)

    if not traits:
        raise AkgwasError(E_INPUT_MISSING, "no sigk files found; run the kgwas/collect stages first")
    aligner = cfg.alignment.aligner
    if aligner not in ("bowtie", "bwa"):
        raise AkgwasError(E_ALIGNER_UNKNOWN, f"unknown aligner: {aligner!r}")

    ensure_dir(os.path.join(final_dir, "samfile"))
    ensure_dir(os.path.join(final_dir, "sigk_mapped"))

    if ctx.dry_run:
        ctx.log(f"DRY-RUN: alignment {len(traits)} trait(s) x {len(refs)} reference(s)")
        return

    require_tools([aligner])
    errors: list[str] = []
    workers = max(1, min(ctx.cores, len(refs)))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_align_reference, ctx, ref, traits): ref for ref in refs}
        for future in as_completed(futures):
            ref = futures[future]
            try:
                future.result()
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{ref.name}: {exc}")
    if errors:
        raise AkgwasError(E_ALIGN_FAILED, "alignment failures:\n  " + "\n  ".join(errors[:10]))

    ctx.markers.mark_done(key, inputs)
    ctx.log("alignment: done")


def _traits(final_dir: str) -> list[str]:
    return sorted(
        os.path.basename(p)[: -len("_sigK.txt")]
        for p in glob.glob(os.path.join(final_dir, "sigk", "*_sigK.txt"))
    )


def _refs(ctx: RunContext) -> list[RefInfo]:
    aligner = ctx.cfg.alignment.aligner
    refs = []
    seen: dict[str, int] = {}
    for entry in read_ref_list(ctx.cfg.ref_list_path):
        index = entry
        if aligner == "bwa" and index.endswith(".sa"):
            index = index[: -len(".sa")]
        name = _strip_seq_ext(os.path.basename(index)) or os.path.basename(index)
        if name in seen:
            seen[name] += 1
            name = f"{name}_{seen[name]}"
        else:
            seen[name] = 0
        refs.append(RefInfo(entry=entry, index=index, name=name))
    return refs


def _strip_seq_ext(name: str) -> str:
    changed = True
    while changed:
        changed = False
        for ext in _SEQ_EXTS:
            if name.endswith(ext):
                name = name[: -len(ext)]
                changed = True
    return name


def _align_reference(ctx: RunContext, ref: RefInfo, traits: list[str]) -> None:
    """All traits against one reference, sequentially (keeps the index hot)."""
    cfg = ctx.cfg
    final_dir = cfg.final_path
    for trait in traits:
        fasta = os.path.join(final_dir, "fasta", f"{trait}_fasta.txt")
        if not os.path.exists(fasta) or os.path.getsize(fasta) == 0:
            ctx.log(f"alignment: {trait} -> {ref.name}: no significant k-mers, skip")
            continue
        sam = os.path.join(final_dir, "samfile", f"{trait}_{ref.name}.sam")
        ctx.run_tool(
            _aligner_cmd(cfg, ref, fasta),
            stdout_path=sam,
            error_code=E_ALIGN_FAILED,
        )
        positions = _parse_sam(sam, cfg.kmer_len, cfg.alignment.max_hits)
        _write_join(ctx, trait, ref, positions)


def _aligner_cmd(cfg, ref: RefInfo, fasta: str) -> list[str]:
    aligner = cfg.alignment.aligner
    threads = cfg.alignment.threads
    if aligner == "bowtie":
        # -a: all; -v 0: end-to-end exact; -f: fasta input; -S: SAM to stdout.
        return ["bowtie", "-a", "-v", "0", "-f", "-S", "-p", threads, ref.index, fasta]
    # bwa mem: -k/-T set to k so only full-length exact seeds pass; -a all; -c cap.
    return [
        "bwa", "mem", "-k", cfg.kmer_len, "-T", cfg.kmer_len, "-a",
        "-c", cfg.alignment.max_hits, "-t", threads, ref.index, fasta,
    ]


def _parse_sam(sam_path: str, k: int, max_hits: int) -> dict[str, list[tuple[str, str, str, str, str]]]:
    """Return {kmer: [(flag, rname, pos, mapq, cigar), ...]} for full-length exact hits."""
    positions: dict[str, list[tuple[str, str, str, str, str]]] = {}
    full_cigar = f"{k}M"
    with open(sam_path, encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            qname, flag_s, rname, pos, mapq, cigar = fields[0], fields[1], fields[2], fields[3], fields[4], fields[5]
            try:
                flag = int(flag_s)
            except ValueError:
                continue
            if flag & 4:  # unmapped
                continue
            if rname == "*":
                continue
            if cigar != full_cigar:
                continue
            nm = _nm_tag(fields[11:])
            if nm is not None and nm != 0:
                continue
            positions.setdefault(qname, []).append((flag_s, rname, pos, mapq, cigar))

    # Drop entire k-mers that exceed the hit cap.
    return {kmer: hits for kmer, hits in positions.items() if len(hits) <= max_hits}


def _nm_tag(tags: list[str]) -> int | None:
    for tag in tags:
        if tag.startswith("NM:i:"):
            try:
                return int(tag[5:])
            except ValueError:
                return None
    return None


def _write_join(ctx: RunContext, trait: str, ref: RefInfo, positions: dict) -> None:
    cfg = ctx.cfg
    sigk = os.path.join(cfg.final_path, "sigk", f"{trait}_sigK.txt")
    values: dict[str, str] = {}
    with open(sigk, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5:
                values[parts[0]] = "\t".join(parts[1:5])

    rows = []
    for kmer, hits in positions.items():
        stat = values.get(kmer)
        if stat is None:
            continue
        for flag, rname, pos, mapq, cigar in hits:
            rows.append((rname, pos, f"{kmer}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{stat}"))
    rows.sort(key=lambda r: (r[0], int(r[1]) if r[1].isdigit() else 0))

    dest = os.path.join(cfg.final_path, "sigk_mapped", f"{trait}_{ref.name}_chrpos_sigk.txt")
    atomic_write_text(dest, "".join(r[2] + "\n" for r in rows))
