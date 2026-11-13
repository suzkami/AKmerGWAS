"""Shared run context passed to every stage."""

from __future__ import annotations

import shlex
from dataclasses import dataclass

from .config import Config
from .logging_utils import Logger
from .tools import MarkerStore, RunResult, run


#: Stages whose results depend on the phenotype (or its derivatives). When a
#: run name is set these get a per-run completion marker so different phenotype
#: runs sharing one outdir do not skip each other. build/multibed/qk are
#: trait-independent and keep a single shared marker.
RUN_SCOPED_STAGES = frozenset({"kgwas", "collect", "alignment"})


@dataclass
class RunContext:
    cfg: Config
    logger: Logger
    markers: MarkerStore
    cores: int = 1
    dry_run: bool = False
    force: bool = False

    def run_tool(self, cmd, *, error_code: str = "E999", **kwargs) -> RunResult:
        """Run an external tool (honouring --dry-run)."""
        cmd = [str(c) for c in cmd]
        if self.dry_run:
            self.logger.info("DRY-RUN: " + " ".join(shlex.quote(c) for c in cmd))
            return RunResult(0)
        return run(cmd, log=self.logger.info, error_code=error_code, **kwargs)

    def log(self, message: str) -> None:
        self.logger.info(message)

    def marker_key(self, stage: str) -> str:
        """Completion-marker key for ``stage`` (run-scoped when a run name is set)."""
        if self.cfg.run_name and stage in RUN_SCOPED_STAGES:
            return f"{self.cfg.run_name}.{stage}"
        return stage
