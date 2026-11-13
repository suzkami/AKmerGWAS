"""Stage ``build``: k-mer counting with kmtricks (count + merge).

Produces the kmtricks run directory. The aggregate (text) step lives in the
``multibed`` stage so its output can be streamed straight into bed chunks
without ever writing the giant text matrix to disk.
"""

from __future__ import annotations

import os

from ..context import RunContext
from ..errors import E_KMTricks_FAILED
from ..tools import ensure_dir, remove_if_exists, require_tools

STAGE = "build"


def run(ctx: RunContext) -> None:
    cfg = ctx.cfg
    run_dir = cfg.kmtricks_run_path
    inputs = [cfg.fof_path]

    if not ctx.force and ctx.markers.is_done(STAGE, inputs):
        ctx.log("build: marker valid -> skip")
        return
    ctx.markers.clear(STAGE)

    if not ctx.dry_run:
        require_tools(["kmtricks"])
        # kmtricks refuses to start when --run-dir already exists, so start from
        # a clean slate (this also makes a --force re-run safe). Create only the
        # parent; kmtricks itself creates the run directory.
        remove_if_exists(run_dir)
        ensure_dir(os.path.dirname(run_dir))

    b = cfg.build
    count_cmd = [
        "kmtricks", "pipeline",
        "--file", cfg.fof_path,
        "--run-dir", run_dir,
        "--kmer-size", cfg.kmer_len,
        "--mode", "kmer:pa:bin",
        "--hard-min", b.hard_min,
        "--recurrence-min", b.recurrence_min,
        "--nb-partitions", b.partitions,
        "-t", b.threads,
        "--until", "count",
        "--cpr",
    ]
    ctx.run_tool(count_cmd, error_code=E_KMTricks_FAILED)

    merge_cmd = [
        "kmtricks", "merge",
        "--run-dir", run_dir,
        "--recurrence-min", b.recurrence_min,
        "--mode", "kmer:pa:bin",
        "--cpr",
        "-t", "2",  # kmtricks merge is the one fixed-thread step (I/O bound)
    ]
    ctx.run_tool(merge_cmd, error_code=E_KMTricks_FAILED)

    if not ctx.dry_run:
        ctx.markers.mark_done(STAGE, inputs)
        ctx.log(f"build: done -> {run_dir}")
