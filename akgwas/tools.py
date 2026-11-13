"""Subprocess wrapper, filesystem helpers, markers and preflight checks.

This module is deliberately dependency-free and is the only place that talks to
the outside world (processes, ``os.replace``, ``resource`` limits). Keeping the
side effects in one module makes the rest of the code easy to reason about.
"""

from __future__ import annotations

import hashlib
import json
import os
import shlex
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

from .errors import (
    AkgwasError,
    E_INTERNAL,
    E_INPUT_MISSING,
    E_TOOL_MISSING,
    E_ULIMIT_TOO_LOW,
    E_WORKDIR_UNWRITABLE,
)

# ---------------------------------------------------------------------------
# Sizes
# ---------------------------------------------------------------------------

_UNITS = {
    "": 1,
    "B": 1,
    "K": 1024,
    "KB": 1024,
    "KIB": 1024,
    "M": 1024 ** 2,
    "MB": 1024 ** 2,
    "MIB": 1024 ** 2,
    "G": 1024 ** 3,
    "GB": 1024 ** 3,
    "GIB": 1024 ** 3,
    "T": 1024 ** 4,
    "TB": 1024 ** 4,
    "TIB": 1024 ** 4,
}


def parse_size(value: str | int) -> int:
    """Parse ``"1G"`` / ``"500M"`` / ``1073741824`` into bytes.

    ``K/M/G/T`` are binary multiples (matching ``parallel --block`` semantics).
    """
    if isinstance(value, int):
        return value
    text = str(value).strip().upper()
    if not text:
        raise AkgwasError(E_INTERNAL, "empty size value")
    num = ""
    unit = ""
    for ch in text:
        if ch.isdigit() or ch == ".":
            num += ch
        else:
            unit += ch
    if not num:
        raise AkgwasError(E_INTERNAL, f"invalid size value: {value!r}")
    if unit not in _UNITS:
        raise AkgwasError(E_INTERNAL, f"unknown size unit in {value!r} (use K/M/G/T)")
    return int(float(num) * _UNITS[unit])


def human_bytes(n: int) -> str:
    for unit, factor in (("T", 1024 ** 4), ("G", 1024 ** 3), ("M", 1024 ** 2), ("K", 1024)):
        if n >= factor:
            return f"{n / factor:.1f}{unit}"
    return f"{n}B"


# ---------------------------------------------------------------------------
# Filesystem
# ---------------------------------------------------------------------------

def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def check_writable(path: str) -> None:
    try:
        ensure_dir(path)
        probe = os.path.join(path, ".akgwas_write_probe")
        with open(probe, "w", encoding="utf-8") as fh:
            fh.write("ok")
        os.remove(probe)
    except OSError as exc:
        raise AkgwasError(E_WORKDIR_UNWRITABLE, f"cannot write to {path}: {exc}") from exc


def atomic_write_text(path: str, text: str) -> None:
    """Write ``text`` to ``path`` atomically (tmp file + ``os.replace``)."""
    ensure_dir(os.path.dirname(os.path.abspath(path)))
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), prefix=".tmp_")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def atomic_replace(tmp_path: str, final_path: str) -> None:
    ensure_dir(os.path.dirname(os.path.abspath(final_path)))
    os.replace(tmp_path, final_path)


def tmp_sibling(path: str) -> str:
    """Return a temp path next to ``path`` (same filesystem, so replace is atomic)."""
    ensure_dir(os.path.dirname(os.path.abspath(path)))
    return path + ".tmp"


def remove_if_exists(path: str) -> None:
    try:
        if os.path.isdir(path):
            shutil.rmtree(path)
        elif os.path.exists(path):
            os.remove(path)
    except OSError:
        pass


# ---------------------------------------------------------------------------
# External tools
# ---------------------------------------------------------------------------

def tool_exists(name: str) -> bool:
    return shutil.which(name) is not None


def require_tools(names: Iterable[str]) -> None:
    missing = [n for n in names if not tool_exists(n)]
    if missing:
        raise AkgwasError(
            E_TOOL_MISSING,
            "required tool(s) not found on PATH: " + ", ".join(missing),
        )


def tool_version(name: str) -> str:
    """Best-effort version string for the run manifest (never raises)."""
    for flag in ("--version", "version", "-v"):
        try:
            proc = subprocess.run(
                [name, flag], capture_output=True, text=True, timeout=10
            )
        except (OSError, subprocess.SubprocessError):
            continue
        out = (proc.stdout or proc.stderr or "").strip().splitlines()
        if out:
            return out[0][:120]
    return "unknown"


@dataclass
class RunResult:
    returncode: int
    stdout: str = ""
    stderr: str = ""


def run(
    cmd: Sequence[str],
    *,
    cwd: str | None = None,
    capture: bool = False,
    check: bool = True,
    env: Mapping[str, str] | None = None,
    stdout_path: str | None = None,
    log=print,
    error_code: str = E_INTERNAL,
) -> RunResult:
    """Run a subprocess, echoing the command.

    ``capture=True`` captures combined output (and still returns it). Otherwise
    output is inherited so long-running tools stream to the console/log.
    """
    cmd = [str(c) for c in cmd]
    log("RUN: " + " ".join(shlex.quote(c) for c in cmd))
    full_env = None
    if env is not None:
        full_env = dict(os.environ)
        full_env.update(env)

    stdout_fh = None
    try:
        if stdout_path is not None:
            stdout_fh = open(stdout_path, "w", encoding="utf-8")
        if capture:
            proc = subprocess.run(
                cmd, cwd=cwd, env=full_env, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            )
            result = RunResult(proc.returncode, proc.stdout or "")
            if result.stdout:
                log(result.stdout.rstrip("\n"))
        else:
            proc = subprocess.run(cmd, cwd=cwd, env=full_env, stdout=stdout_fh)
            result = RunResult(proc.returncode)
    except FileNotFoundError as exc:
        raise AkgwasError(E_TOOL_MISSING, f"command not found: {cmd[0]}") from exc
    finally:
        if stdout_fh is not None:
            stdout_fh.close()

    if check and result.returncode != 0:
        raise AkgwasError(
            error_code,
            f"command failed (exit {result.returncode}): " + " ".join(shlex.quote(c) for c in cmd),
        )
    return result


def popen_stdout(cmd: Sequence[str], *, cwd: str | None = None, env: Mapping[str, str] | None = None):
    """Start a process and return it with ``stdout`` as a text pipe.

    Used to stream large kmtricks text output without ever writing it to disk.
    ``stderr`` is redirected to an unnamed temporary file (not a pipe): the
    caller reads the pipe only after the process exits, so a pipe there could
    fill up and deadlock a tool that logs progress to stderr. The file is
    attached to the returned ``Popen`` as ``akgwas_stderr``; read it with
    :func:`popen_stderr`.
    """
    cmd = [str(c) for c in cmd]
    full_env = None
    if env is not None:
        full_env = dict(os.environ)
        full_env.update(env)
    err = tempfile.TemporaryFile("w+", encoding="utf-8")
    proc = subprocess.Popen(
        cmd, cwd=cwd, env=full_env, text=True,
        stdout=subprocess.PIPE, stderr=err,
    )
    proc.akgwas_stderr = err  # type: ignore[attr-defined]
    return proc


def popen_stderr(proc) -> str:
    """Read and close the temporary stderr file attached by :func:`popen_stdout`."""
    err = getattr(proc, "akgwas_stderr", None)
    if err is None:
        return ""
    try:
        err.seek(0)
        return err.read()
    finally:
        err.close()
        proc.akgwas_stderr = None


# ---------------------------------------------------------------------------
# ulimit
# ---------------------------------------------------------------------------

def current_nofile_limit() -> int | None:
    try:
        import resource  # noqa: PLC0415  (Unix only; missing on Windows)

        soft, _hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        return None if soft == resource.RLIM_INFINITY else int(soft)
    except (ImportError, ValueError, OSError):
        return None


def check_ulimit_n(minimum: int) -> None:
    current = current_nofile_limit()
    if current is not None and current < minimum:
        raise AkgwasError(
            E_ULIMIT_TOO_LOW,
            f"open-file limit is {current}, need >= {minimum}",
        )


# ---------------------------------------------------------------------------
# Completion markers + fingerprints
# ---------------------------------------------------------------------------

def fingerprint(paths: Iterable[str]) -> str:
    """Hash of (abspath, size, mtime_ns) for each input path.

    Directories are expanded one level so a chunk directory changes when files
    are added/removed. Missing paths contribute a ``MISSING`` token so a
    deleted input invalidates the marker.
    """
    items: list[tuple[str, str, str]] = []
    for raw in paths:
        if raw is None:
            continue
        path = os.path.abspath(raw)
        if os.path.isdir(path):
            for entry in sorted(os.listdir(path)):
                fp = os.path.join(path, entry)
                try:
                    st = os.stat(fp)
                    items.append((fp, str(st.st_size), str(st.st_mtime_ns)))
                except OSError:
                    items.append((fp, "MISSING", "0"))
        else:
            try:
                st = os.stat(path)
                items.append((path, str(st.st_size), str(st.st_mtime_ns)))
            except OSError:
                items.append((path, "MISSING", "0"))
    blob = json.dumps(sorted(items), separators=(",", ":"))
    return hashlib.sha1(blob.encode("utf-8")).hexdigest()


class MarkerStore:
    """Tracks per-stage completion using atomic ``.done`` files under ``.akgwas``."""

    def __init__(self, workdir: str, version: str) -> None:
        self.dir = os.path.join(workdir, ".akgwas")
        self.version = version

    def _path(self, stage: str) -> str:
        return os.path.join(self.dir, f"{stage}.done")

    def is_done(self, stage: str, inputs: Iterable[str]) -> bool:
        path = self._path(stage)
        if not os.path.exists(path):
            return False
        try:
            with open(path, encoding="utf-8") as fh:
                data = json.load(fh)
        except (OSError, json.JSONDecodeError):
            return False
        if data.get("version") != self.version:
            return False
        return data.get("fingerprint") == fingerprint(inputs)

    def mark_done(self, stage: str, inputs: Iterable[str], extra: Mapping | None = None) -> None:
        payload = {
            "stage": stage,
            "version": self.version,
            "fingerprint": fingerprint(inputs),
            "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
        }
        if extra:
            payload.update(dict(extra))
        atomic_write_text(self._path(stage), json.dumps(payload, indent=2))

    def clear(self, stage: str) -> None:
        remove_if_exists(self._path(stage))


# ---------------------------------------------------------------------------
# Small IO helpers
# ---------------------------------------------------------------------------

def read_lines(path: str) -> list[str]:
    with open(path, encoding="utf-8") as fh:
        return [ln.rstrip("\n") for ln in fh]


def require_files(paths: Iterable[str]) -> None:
    for path in paths:
        if not os.path.exists(path):
            raise AkgwasError(E_INPUT_MISSING, f"required file not found: {path}")
