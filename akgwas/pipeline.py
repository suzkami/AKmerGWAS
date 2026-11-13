"""Pipeline orchestration: preflight, fixed DAG, markers, manifest.

Replaces the old ``cli.py -> *.smk -> *.sh`` stack with a single Python driver.
The DAG is linear and fixed, so no scheduler is needed; the only fan-out points
(``chunk x trait`` in kgwas, ``reference`` in alignment) open their own thread
pools inside the stages.
"""

from __future__ import annotations

import json
import os
import time

from . import __version__
from .config import Config
from .context import RunContext
from .errors import AkgwasError, E_CONFIG_INVALID, E_ROW_COUNT_MISMATCH, E_SAMPLE_ID_MISMATCH
from .inputs import read_fof_samples, read_ref_list
from .logging_utils import Logger
from .pheno import read_pheno
from .stages import alignment, build, collect, kgwas, multibed, qk
from .tools import (
    MarkerStore,
    check_ulimit_n,
    check_writable,
    ensure_dir,
    require_files,
    require_tools,
    tool_exists,
    tool_version,
)

# module name -> ordered list of stages
MODULE_STAGES: dict[str, list[str]] = {
    "build": ["build", "multibed"],
    "kgwas": ["qk", "kgwas", "collect"],
    "alignment": ["alignment"],
    "all": ["build", "multibed", "qk", "kgwas", "collect", "alignment"],
}

_STAGE_FUNCS = {
    "build": build.run,
    "multibed": multibed.run,
    "qk": qk.run,
    "kgwas": kgwas.run,
    "collect": collect.run,
    "alignment": alignment.run,
}


def run_module(
    cfg: Config,
    module: str,
    *,
    cores: int = 1,
    dry_run: bool = False,
    force: bool = False,
    verbose: bool = True,
) -> int:
    if module not in MODULE_STAGES:
        raise AkgwasError(E_CONFIG_INVALID, f"unknown module: {module!r}")

    stages = MODULE_STAGES[module]

    ensure_dir(cfg.outdir)
    check_writable(cfg.outdir)
    markers = MarkerStore(cfg.outdir, __version__)
    logger = Logger(os.path.join(cfg.outdir, ".akgwas", "akgwas.log"), verbose=verbose)
    ctx = RunContext(cfg=cfg, logger=logger, markers=markers, cores=cores, dry_run=dry_run, force=force)

    started = time.time()
    logger.phase(f"akgwas {__version__} :: {module} :: {', '.join(stages)}")
    if cfg.run_name:
        logger.info(f"run name: {cfg.run_name} -> {cfg.run_root}")
    try:
        preflight(ctx, stages)
        _write_manifest(ctx, module)
        for stage in stages:
            _STAGE_FUNCS[stage](ctx)
        logger.info(f"module {module} finished in {time.time() - started:.1f}s")
        return 0
    except AkgwasError as exc:
        # Always surfaced (even under --quiet): logger.error() echoes to stderr
        # unconditionally, so a failure is never silent.
        logger.error(exc.render())
        return 1
    finally:
        logger.close()


def preflight(ctx: RunContext, stages: list[str]) -> None:
    """Fail fast, with a clear message, before touching any data."""
    cfg = ctx.cfg
    if not ctx.dry_run:
        require_tools(_stage_tools(cfg, stages))
        check_ulimit_n(cfg.ulimit_min)

    require_files([cfg.fof_path])

    # Inputs relevant to the requested module.
    if any(s in stages for s in ("kgwas", "collect")):
        require_files([cfg.pheno_path])
        _check_sample_consistency(cfg)
    if "alignment" in stages:
        require_files([cfg.ref_list_path])
        read_ref_list(cfg.ref_list_path)


def _stage_tools(cfg: Config, stages: list[str]) -> list[str]:
    tools: list[str] = []
    if "build" in stages or "multibed" in stages:
        tools.append("kmtricks")
    if "multibed" in stages and cfg.multibed.backend == "plink":
        tools.append("plink")
    if "qk" in stages and not cfg.qk.use_external:
        tools += ["plink", "gemma"]
    if "kgwas" in stages:
        tools.append("gemma")
    if "alignment" in stages:
        tools.append(cfg.alignment.aligner)
    # de-dup, preserve order
    seen = set()
    return [t for t in tools if not (t in seen or seen.add(t))]


def _check_sample_consistency(cfg: Config) -> None:
    """Early cross-check of fof sample names vs phenotype sample IDs."""
    fof_ids = read_fof_samples(cfg.fof_path)
    table = read_pheno(cfg.pheno_path)
    if len(table.ids) != len(fof_ids):
        raise AkgwasError(
            E_ROW_COUNT_MISMATCH,
            f"fof has {len(fof_ids)} samples but phenotype has {len(table.ids)} rows",
        )
    if set(fof_ids) != set(table.ids):
        missing = sorted(set(fof_ids) - set(table.ids))[:5]
        extra = sorted(set(table.ids) - set(fof_ids))[:5]
        raise AkgwasError(
            E_SAMPLE_ID_MISMATCH,
            f"fof and phenotype sample IDs differ (missing={missing}, extra={extra})",
        )


def _write_manifest(ctx: RunContext, module: str) -> None:
    cfg = ctx.cfg
    if ctx.dry_run:
        return
    tool_names = ["kmtricks", "gemma", "plink", "bwa", "bowtie"]
    versions = {t: tool_version(t) for t in tool_names if tool_exists(t)}
    manifest = {
        "akgwas_version": __version__,
        "module": module,
        "cores": ctx.cores,
        "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "config": cfg.to_manifest(),
        "tools": versions,
    }
    ensure_dir(os.path.dirname(cfg.manifest_path))
    with open(cfg.manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
