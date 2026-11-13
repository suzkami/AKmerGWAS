"""A tiny dependency-free logger.

AKmerGWAS runs many external tools; the logger gives each run a single log file
(``<workdir>/.akgwas/akgwas.log``) plus console output, with timestamps and
phase/step prefixes so a long multi-hour run is still readable afterwards.

Errors are special-cased: :meth:`Logger.error` always echoes to stderr, even
under ``--quiet`` / ``verbose=False``. A silent failure is the hardest kind to
debug, and the pipeline reports user-facing failures through this one method.
"""

from __future__ import annotations

import os
import sys
import time
from typing import TextIO


class Logger:
    def __init__(self, logfile: str | None = None, *, verbose: bool = True) -> None:
        self.verbose = verbose
        self._fh: TextIO | None = None
        self._start = time.time()
        if logfile:
            os.makedirs(os.path.dirname(os.path.abspath(logfile)), exist_ok=True)
            self._fh = open(logfile, "a", encoding="utf-8", buffering=1)

    def _emit(self, level: str, message: str, *, console: bool = False) -> None:
        elapsed = time.time() - self._start
        line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')} +{elapsed:8.1f}s] {level:<5} {message}"
        if self.verbose or console:
            print(line, file=sys.stderr, flush=True)
        if self._fh is not None:
            self._fh.write(line + "\n")

    def debug(self, message: str) -> None:
        self._emit("DEBUG", message)

    def info(self, message: str) -> None:
        self._emit("INFO", message)

    def warn(self, message: str) -> None:
        self._emit("WARN", message)

    def error(self, message: str) -> None:
        # Always shown (not gated by --verbose): failures must never be silent.
        self._emit("ERROR", message, console=True)

    def phase(self, name: str) -> None:
        self._emit("PHASE", f"=== {name} ===")

    def close(self) -> None:
        if self._fh is not None:
            self._fh.close()
            self._fh = None
