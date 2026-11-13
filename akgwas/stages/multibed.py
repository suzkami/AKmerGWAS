"""Stage ``multibed``: stream the kmtricks text matrix into PLINK bed chunks.

This replaces the old ``GNU parallel --block`` + ``awk`` + ``plink`` pipeline.
Key properties:

* the kmtricks aggregate text is consumed as a **stream** (stdout pipe, or a
  temporary file that is deleted) -- the multi-gigabyte ``Origin_matrix.txt``
  and all ``.tped`` intermediates are gone;
* chunks are cut by a **byte budget** (``block_size``), guarded by
  ``min_block_size`` / ``max_block_size``, so the chunk count is controlled;
* two backends: ``python`` writes ``.bed`` directly; ``plink`` writes one
  temporary ``.tped`` per chunk and lets plink convert (deleted immediately, so
  peak disk stays ~one chunk).

MAF filter (kept identical to the reviewed-correct original): a k-mer survives
iff ``max(#absent, #present) <= int(n - maf*n)``.
"""

from __future__ import annotations

import glob
import os
from typing import Iterator, Sequence

from ..bed import BedWriter, HOM_A1, HOM_A2, bytes_per_variant, pack_variant
from ..context import RunContext
from ..errors import (
    AkgwasError,
    E_BED_BACKEND_UNKNOWN,
    E_EMPTY_MATRIX,
    E_KMTricks_FAILED,
    E_MATRIX_MALFORMED,
)
from ..inputs import read_fof_samples
from ..tools import (
    atomic_write_text,
    ensure_dir,
    parse_size,
    popen_stderr,
    popen_stdout,
    remove_if_exists,
    require_tools,
)

STAGE = "multibed"


def run(ctx: RunContext) -> None:
    cfg = ctx.cfg
    run_dir = cfg.kmtricks_run_path
    bed_dir = cfg.bed_path
    inputs = [run_dir, cfg.fof_path]

    if not ctx.force and ctx.markers.is_done(STAGE, inputs):
        ctx.log("multibed: marker valid -> skip")
        return
    ctx.markers.clear(STAGE)

    if cfg.multibed.backend == "plink":
        require_tools(["kmtricks", "plink"])
    elif cfg.multibed.backend == "python":
        require_tools(["kmtricks"])
    else:
        raise AkgwasError(
            E_BED_BACKEND_UNKNOWN, f"unknown bed backend: {cfg.multibed.backend!r}"
        )

    ensure_dir(bed_dir)
    # Clean any stale chunks so a re-run never mixes old and new output.
    for stale in glob.glob(os.path.join(bed_dir, "matrix_*")):
        remove_if_exists(stale)

    sample_ids = read_fof_samples(cfg.fof_path)
    n = len(sample_ids)
    fam_lines = [f"0 {sid} 0 0 0 0" for sid in sample_ids]
    atomic_write_text(
        os.path.join(bed_dir, "sample_id_inorder.txt"), "\n".join(sample_ids) + "\n"
    )

    if ctx.dry_run:
        ctx.log(f"DRY-RUN: multibed would stream kmtricks aggregate -> {bed_dir}/matrix_*.bed")
        return

    stats = _stream_to_chunks(ctx, _matrix_line_iter(ctx), bed_dir, fam_lines, n)
    if stats["variants"] == 0:
        raise AkgwasError(
            E_EMPTY_MATRIX,
            "no k-mers survived the MAF filter; check --hard-min / --recurrence-min and --maf",
        )
    ctx.log(
        f"multibed: {stats['variants']} k-mers in {stats['chunks']} chunk(s); "
        f"{stats['skipped']} filtered by MAF"
    )
    ctx.markers.mark_done(STAGE, inputs, extra={"chunks": stats["chunks"]})


# ---------------------------------------------------------------------------
# Streaming the kmtricks text matrix
# ---------------------------------------------------------------------------

def _matrix_line_iter(ctx: RunContext) -> Iterator[str]:
    cfg = ctx.cfg
    base = [
        "kmtricks", "aggregate",
        "--run-dir", cfg.kmtricks_run_path,
        "--pa-matrix", "kmer",
        "--cpr-in",
        "-t", cfg.multibed.threads,
        "--format", "text",
    ]
    if cfg.build.aggregate_to_stdout:
        proc = popen_stdout(base)
        try:
            assert proc.stdout is not None
            for line in proc.stdout:
                yield line
            rc = proc.wait()
            if rc != 0:
                err = popen_stderr(proc)
                raise AkgwasError(
                    E_KMTricks_FAILED,
                    f"kmtricks aggregate (stdout) exited {rc}: {err.strip()[:300]}",
                )
        finally:
            if proc.poll() is None:
                proc.kill()
            popen_stderr(proc)  # ensure the stderr temp file is closed
    else:
        tmp = os.path.join(cfg.kmtricks_run_path, ".aggregate_text.tmp")
        ctx.run_tool(base + ["--output", tmp], error_code=E_KMTricks_FAILED)
        try:
            with open(tmp, encoding="utf-8") as fh:
                yield from fh
        finally:
            remove_if_exists(tmp)


def _stream_to_chunks(
    ctx: RunContext,
    lines: Iterator[str],
    bed_dir: str,
    fam_lines: Sequence[str],
    n: int,
) -> dict[str, int]:
    cfg = ctx.cfg
    mb = cfg.multibed
    target = parse_size(mb.block_size)
    stride = bytes_per_variant(n)
    limit = int(n - mb.maf * n)

    index = 0
    variants = 0
    skipped = 0
    first_line = True

    sink = None
    sink_bytes = 0

    def open_sink():
        nonlocal index
        index += 1
        prefix = os.path.join(bed_dir, f"matrix_{index:04d}")
        return _make_sink(mb.backend, prefix, fam_lines, n, ctx), prefix

    try:
        for raw in lines:
            fields = raw.split()
            if not fields:
                continue
            if len(fields) != n + 1:
                if first_line:
                    first_line = False
                    ctx.log("multibed: skipping header line from aggregate output")
                    continue
                raise AkgwasError(
                    E_MATRIX_MALFORMED,
                    f"matrix row has {len(fields)} fields, expected {n + 1} (k-mer + {n} samples)",
                )
            first_line = False
            kmer, vals = fields[0], fields[1:]
            zeros = vals.count("0")
            ones = vals.count("1")
            if zeros + ones != n:
                raise AkgwasError(
                    E_MATRIX_MALFORMED,
                    f"matrix row {kmer!r} contains values other than 0/1",
                )
            if not (zeros <= limit and ones <= limit):
                skipped += 1
                continue

            if sink is None or sink_bytes >= target:
                if sink is not None:
                    sink.close()
                sink, _prefix = open_sink()
                sink_bytes = 0
            sink.add(kmer, vals)
            sink_bytes += stride + len(kmer) + 8
            variants += 1
    finally:
        if sink is not None:
            sink.close()

    return {"variants": variants, "skipped": skipped, "chunks": index}


# ---------------------------------------------------------------------------
# Chunk sinks
# ---------------------------------------------------------------------------

def _minor_first_alleles(vals: Sequence[str], n: int) -> tuple[str, str]:
    """Mimic plink's default (minor allele first) so GEMMA sees the same coding."""
    present = vals.count("1")
    absent = n - present
    minor = "1" if absent <= present else "2"
    major = "2" if minor == "1" else "1"
    return minor, major


class _PythonSink:
    """Writes a chunk's ``.bed``/``.bim``/``.fam`` directly (no ``.tped``)."""

    def __init__(self, prefix: str, fam_lines: Sequence[str], n: int, ctx: RunContext) -> None:
        self.n = n
        self.writer = BedWriter(prefix, fam_lines)

    def add(self, kmer: str, vals: Sequence[str]) -> None:
        a1, a2 = _minor_first_alleles(vals, self.n)
        codes = []
        for v in vals:
            is_a1 = (v == "0" and a1 == "1") or (v == "1" and a1 == "2")
            codes.append(HOM_A1 if is_a1 else HOM_A2)
        self.writer.add_variant(["1", kmer, "0", "0", a1, a2], pack_variant(codes))

    def close(self) -> None:
        self.writer.close()


class _PlinkSink:
    """Writes one temporary ``.tped``/``.tfam`` and converts it with plink."""

    def __init__(self, prefix: str, fam_lines: Sequence[str], n: int, ctx: RunContext) -> None:
        self.prefix = prefix
        self.ctx = ctx
        self.n = n
        self.tped_prefix = prefix + ".tpedtmp"
        self._tped = open(self.tped_prefix + ".tped", "w", encoding="utf-8")
        with open(self.tped_prefix + ".tfam", "w", encoding="utf-8") as fh:
            for line in fam_lines:
                fh.write(line + "\n")

    def add(self, kmer: str, vals: Sequence[str]) -> None:
        geno = " ".join("2 2" if v == "1" else "1 1" for v in vals)
        self._tped.write(f"1\t{kmer}\t0\t0\t{geno}\n")

    def close(self) -> None:
        self._tped.close()
        try:
            self.ctx.run_tool(
                ["plink", "--tfile", self.tped_prefix, "--make-bed", "--silent",
                 "--out", self.prefix],
                error_code=E_KMTricks_FAILED,
            )
        finally:
            for ext in ("tped", "tfam", "nosex", "log"):
                remove_if_exists(f"{self.tped_prefix}.{ext}")
            for ext in ("nosex", "log"):
                remove_if_exists(f"{self.prefix}.{ext}")


def _make_sink(backend: str, prefix: str, fam_lines, n, ctx):
    if backend == "python":
        return _PythonSink(prefix, fam_lines, n, ctx)
    if backend == "plink":
        return _PlinkSink(prefix, fam_lines, n, ctx)
    raise AkgwasError(E_BED_BACKEND_UNKNOWN, f"unknown bed backend: {backend!r}")
