"""Run configuration: defaults, CLI overrides and path resolution.

AKmerGWAS is driven **entirely from the command line** -- there is no config
file. Every value has a default (see the dataclasses below); pass the matching
flag to override it. ``programhome`` is the only required value and is the base
for resolving relative paths, so a run is reproducible from the command line
alone (the exact invocation is recorded in the run manifest).

Path rule: ``programhome`` must be absolute; every other path-like value is
resolved relative to it unless already absolute. This removes the old
"relative paths depend on the current working directory" bug -- the same
command works from anywhere.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any

from .errors import AkgwasError, E_CONFIG_INVALID
from .tools import parse_size


@dataclass
class BuildConfig:
    hard_min: int = 2
    recurrence_min: int = 2
    threads: int = 8
    partitions: int = 0
    aggregate_to_stdout: bool = True


@dataclass
class MultibedConfig:
    block_size: str = "1G"
    min_block_size: str = "16M"
    max_block_size: str = "10G"
    maf: float = 0.01
    threads: int = 8
    backend: str = "python"  # "python" (direct bed) | "plink" (tped -> plink)


@dataclass
class QKConfig:
    use_external: bool = False
    kinship: str = ""
    pca: str = ""
    sample_fraction: float = 0.2
    seed: int = 123456789
    gk: int = 1
    pca_has_id: bool = False
    threads: int = 1


@dataclass
class KGwasConfig:
    # One GEMMA call per (chunk, trait) is single-threaded by design; concurrency
    # comes from --cores, so there is no per-job thread knob here.
    threshold: float = 0.0001


@dataclass
class AlignmentConfig:
    ref_list: str = "ref_list.txt"
    aligner: str = "bowtie"  # bowtie | bwa
    max_hits: int = 5000
    threads: int = 1


@dataclass
class CollectConfig:
    # The per-trait, all-k-mers association dump (``gemma_combine/``) is the
    # largest output of a run. The significant k-mers -- with their af/beta/se/
    # p_wald -- are already in ``sigk/`` and ``sigk_mapped/``, so the full dump
    # is only needed to re-threshold or plot, and is therefore opt-in.
    gemma_combine: bool = False


@dataclass
class Config:
    # [common]
    programhome: str = ""
    outdir: str = ""            # defaults to programhome (see make_config)
    run_name: str = ""          # optional: isolates per-phenotype GWAS outputs
    fof_file: str = "fof.txt"
    pheno_file: str = "pheno.txt"
    kmer_len: int = 31
    ulimit_min: int = 4096
    # output directory names (relative to outdir); the first group is shared
    # across runs, the second group is per-run (under run_root).
    bed_dir: str = "bed"
    qk_dir: str = "QK"
    temp_dir: str = "temp_gemma"
    final_dir: str = "kmergwas_outputs"
    kmtricks_run_dir: str = "kmatrixs_out"

    build: BuildConfig = field(default_factory=BuildConfig)
    multibed: MultibedConfig = field(default_factory=MultibedConfig)
    qk: QKConfig = field(default_factory=QKConfig)
    kgwas: KGwasConfig = field(default_factory=KGwasConfig)
    alignment: AlignmentConfig = field(default_factory=AlignmentConfig)
    collect: CollectConfig = field(default_factory=CollectConfig)

    # -- resolved absolute paths ------------------------------------------
    def _resolve(self, value: str) -> str:
        if not value:
            return ""
        if os.path.isabs(value):
            return os.path.normpath(value)
        return os.path.normpath(os.path.join(self.programhome, value))

    @property
    def fof_path(self) -> str:
        return self._resolve(self.fof_file)

    @property
    def pheno_path(self) -> str:
        return self._resolve(self.pheno_file)

    @property
    def ref_list_path(self) -> str:
        return self._resolve(self.alignment.ref_list)

    @property
    def external_kinship_path(self) -> str:
        return self._resolve(self.qk.kinship)

    @property
    def external_pca_path(self) -> str:
        return self._resolve(self.qk.pca)

    # Shared (trait-independent) outputs live under outdir.
    @property
    def bed_path(self) -> str:
        return os.path.join(self.outdir, self.bed_dir)

    @property
    def qk_path(self) -> str:
        return os.path.join(self.outdir, self.qk_dir)

    @property
    def kmtricks_run_path(self) -> str:
        return os.path.join(self.outdir, self.kmtricks_run_dir)

    # Per-run outputs live under run_root (outdir, or outdir/<run_name>).
    @property
    def run_root(self) -> str:
        if self.run_name:
            return os.path.join(self.outdir, self.run_name)
        return self.outdir

    @property
    def pheno_single_path(self) -> str:
        return os.path.join(self.run_root, "pheno_single")

    @property
    def temp_path(self) -> str:
        return os.path.join(self.run_root, self.temp_dir)

    @property
    def final_path(self) -> str:
        return os.path.join(self.run_root, self.final_dir)

    @property
    def manifest_path(self) -> str:
        stem = f"run_manifest_{self.run_name}.json" if self.run_name else "run_manifest.json"
        return os.path.join(self.outdir, ".akgwas", stem)

    def to_manifest(self) -> dict[str, Any]:
        return {
            "programhome": self.programhome,
            "outdir": self.outdir,
            "run_name": self.run_name,
            "fof_file": self.fof_path,
            "pheno_file": self.pheno_path,
            "kmer_len": self.kmer_len,
            "ulimit_min": self.ulimit_min,
            "paths": {
                "bed_dir": self.bed_dir,
                "qk_dir": self.qk_dir,
                "temp_dir": self.temp_dir,
                "final_dir": self.final_dir,
                "kmtricks_run_dir": self.kmtricks_run_dir,
            },
            "build": vars(self.build),
            "multibed": vars(self.multibed),
            "qk": vars(self.qk),
            "kgwas": vars(self.kgwas),
            "alignment": vars(self.alignment),
            "collect": vars(self.collect),
        }


def make_config(programhome: str) -> Config:
    """A fresh config rooted at an absolute ``programhome`` (outdir defaults there)."""
    cfg = Config()
    cfg.programhome = os.path.normpath(os.path.abspath(programhome))
    cfg.outdir = cfg.programhome
    return cfg


#: CLI dest -> (section, attribute). Only *set* destinations override.
OVERRIDES = {
    "fof": ("common", "fof_file"),
    "pheno": ("common", "pheno_file"),
    "kmer_len": ("common", "kmer_len"),
    "ulimit_min": ("common", "ulimit_min"),
    "bed_dir": ("common", "bed_dir"),
    "qk_dir": ("common", "qk_dir"),
    "temp_dir": ("common", "temp_dir"),
    "final_dir": ("common", "final_dir"),
    "kmtricks_run_dir": ("common", "kmtricks_run_dir"),
    "hard_min": ("build", "hard_min"),
    "recurrence_min": ("build", "recurrence_min"),
    "build_threads": ("build", "threads"),
    "partitions": ("build", "partitions"),
    "aggregate_to_stdout": ("build", "aggregate_to_stdout"),
    "block_size": ("multibed", "block_size"),
    "min_block_size": ("multibed", "min_block_size"),
    "max_block_size": ("multibed", "max_block_size"),
    "maf": ("multibed", "maf"),
    "multibed_threads": ("multibed", "threads"),
    "bed_backend": ("multibed", "backend"),
    "use_external_qk": ("qk", "use_external"),
    "kinship": ("qk", "kinship"),
    "pca": ("qk", "pca"),
    "sample_fraction": ("qk", "sample_fraction"),
    "seed": ("qk", "seed"),
    "gk": ("qk", "gk"),
    "pca_has_id": ("qk", "pca_has_id"),
    "qk_threads": ("qk", "threads"),
    "threshold": ("kgwas", "threshold"),
    "ref_list": ("alignment", "ref_list"),
    "aligner": ("alignment", "aligner"),
    "max_hits": ("alignment", "max_hits"),
    "alignment_threads": ("alignment", "threads"),
    "gemma_combine": ("collect", "gemma_combine"),
}


def apply_overrides(cfg: Config, args) -> None:
    """Copy any *explicitly provided* CLI values from ``args`` onto ``cfg``."""
    if getattr(args, "outdir", None):
        cfg.outdir = cfg._resolve(args.outdir)
    if getattr(args, "run_name", None):
        cfg.run_name = args.run_name
    for dest, (section, attr) in OVERRIDES.items():
        value = getattr(args, dest, None)
        if value is None:
            continue
        if section == "common":
            setattr(cfg, attr, value)
        else:
            setattr(getattr(cfg, section), attr, value)


def config_from_args(args) -> Config:
    """Build a validated :class:`Config` from parsed CLI arguments."""
    programhome = getattr(args, "programhome", None)
    if not programhome:
        raise AkgwasError(E_CONFIG_INVALID, "--programhome is required")
    cfg = make_config(programhome)
    apply_overrides(cfg, args)
    validate(cfg)
    return cfg


def validate(cfg: Config) -> None:
    if not cfg.programhome:
        raise AkgwasError(E_CONFIG_INVALID, "programhome is required")
    if not os.path.isabs(cfg.programhome):
        raise AkgwasError(E_CONFIG_INVALID, "programhome must be an absolute path")
    if cfg.kmer_len <= 0:
        raise AkgwasError(E_CONFIG_INVALID, "kmer_len must be positive")
    if cfg.multibed.backend not in ("python", "plink"):
        raise AkgwasError(E_CONFIG_INVALID, "bed backend must be 'python' or 'plink'")
    if cfg.alignment.aligner not in ("bowtie", "bwa"):
        raise AkgwasError(E_CONFIG_INVALID, "aligner must be 'bowtie' or 'bwa'")
    if not 0.0 <= cfg.multibed.maf < 0.5:
        raise AkgwasError(E_CONFIG_INVALID, "maf must be in [0, 0.5)")
    if not 0.0 < cfg.qk.sample_fraction <= 1.0:
        raise AkgwasError(E_CONFIG_INVALID, "sample_fraction must be in (0, 1]")
    lo = parse_size(cfg.multibed.min_block_size)
    hi = parse_size(cfg.multibed.max_block_size)
    target = parse_size(cfg.multibed.block_size)
    if not lo <= target <= hi:
        raise AkgwasError(
            E_CONFIG_INVALID,
            f"block_size ({cfg.multibed.block_size}) must lie within "
            f"[min_block_size, max_block_size] = [{cfg.multibed.min_block_size}, "
            f"{cfg.multibed.max_block_size}]",
        )
