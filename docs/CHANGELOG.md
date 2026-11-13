# Changelog

All notable changes to AKmerGWAS are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

## [1.0.0] — 2026-10-12

Ground-up rewrite of the 0.1.0 Python/Snakemake/Bash stack into a single
pure-Python orchestrator with **no third-party runtime dependencies**
(standard library only). External tools (kmtricks, gemma, plink, bowtie, bwa)
are still invoked as subprocesses.

### Added

- Four subcommands: `akgwas build`, `akgwas kgwas`, `akgwas alignment`,
  `akgwas all`, plus `--dry-run`, `--force`, `--version`, `--errors`.
- A flag for every configuration value, each with a default; only
  `--programhome` is required. `akgwas <subcommand> --help` lists them all,
  grouped by area.
- `--run-name` to isolate phenotype-dependent outputs under `outdir/<name>/`,
  so one `bed/`+`QK/` (computed once) can serve several phenotype runs without
  clobbering or wrongly skipping stages.
- Centralised, stable error codes (`E1xx`–`E9xx`) with human-readable hints and a
  `docs/troubleshooting.md` reference.
- Atomic completion markers (`<stage>.done` + input fingerprint) so re-runs skip
  unchanged stages; `--force` to bypass. Trait-dependent stages are marked
  per run name.
- A run manifest (`.akgwas/run_manifest[_<run>].json`) recording the tool
  version, the full resolved configuration, and external tool versions.
- Two multibed backends: `python` (writes `.bed` directly, default) and `plink`
  (per-chunk `.tped` → plink).
- Selectable aligner: bowtie (`-a -v 0`) or bwa (`mem -k31 -T31 -a -c`).
- Deterministic chunk sampling for the internal K/PCA (`--seed`).
- `--gemma-combine` (default **off**): the per-trait, all-k-mers association
  dump (`kmergwas_outputs/gemma_combine/`) is now opt-in. It is the largest
  output of a run, while the significant k-mers -- with their `af`/`beta`/`se`/
  `p_wald` -- are already in `sigk/` and `sigk_mapped/`, so the full dump is only
  needed to re-threshold or plot. When off, the `collect` stage also skips
  *reading* every chunk's `.assoc.txt` (a large I/O cost in its own right).

### Changed

- **Config: `config.yaml` → command-line flags (no config file).** Every value
  is now a flag with a default; the required `--programhome` is the base for
  resolving relative paths and where outputs are written. The exact invocation
  is recorded in the run manifest. This makes an interim `config.toml` approach
  unnecessary, so v1.0.0 ships with **no configuration file at all**.
- **Phenotype: one single-column file per trait → one multi-column table.** The
  first column is the sample ID, the header row is the trait names, missing
  values are `NA`. The tool validates the ID set against the bed `.fam` and
  reorders rows automatically, then exports the per-trait single-column files
  GEMMA needs. Removes the old `pheno_dir/` layout and its off-by-one bugs.
- **Orchestration: Snakemake `.smk` + GNU parallel/rush + Bash → pure Python.**
  Drops the `snakemake`, `perl`, `parallel`, and `rush` dependencies.
- **Aligner default: `bwa mem` → `bowtie -a -v 0`** (all full-length exact
  hits). bwa remains available.
- **CLI: `--build/--kgwas/--all` flags → subcommands** (`akgwas build|…`).
- **Output layout: shared vs per-run.** `bed/`, `QK/`, `kmatrixs_out/` are
  trait-independent and shared; `pheno_single/`, `temp_gemma/`,
  `kmergwas_outputs/` and their markers are per-run (under `--run-name` when set).
- Sample alignment is now enforced, not advisory: internal K/PCA is built by
  concatenating bed chunks (never `plink --merge-list`, which reorders samples)
  and the PCA sample IDs are checked against the canonical `.fam` order.
- Reference index list: with the bowtie default, `ref_list.txt` holds bowtie
  index **prefixes** (from `bowtie-build`), not `.sa` paths.
- **README Quick start now shows the `bowtie-build` step.** It previously jumped
  straight to `akgwas all` with only a note that `ref_list.txt` holds bowtie
  index prefixes, so a first-time user had no hint that the index has to be
  built first. The Quick start (and the `ref_list.txt` section) now show the
  `bowtie-build <ref.fasta> <prefix>` command and spell out that akgwas never
  builds the index itself — you prepare it and write its prefix into
  `ref_list.txt`.

### Fixed

- **`qk` could save GEMMA's *log* as the kinship matrix.** `_collect_gemma_outputs`
  globbed `<prefix>.*.txt` across both `qk_dir/output/` and `qk_dir`, then took
  the lexically-first hit. On a second `qk` run the previous run's leftover
  `QK_kinship.log.txt` (in `qk_dir`, and 'Q' < 'o') sorted ahead of the new
  `output/QK_kinship.cXX.txt` and was moved over `QK_kinship.txt`. GEMMA then
  rejected that file downstream ("number of columns in the kinship file does not
  match the number of individuals for row = 0") behind a misleading
  `E602`/`.assoc.txt` hint. The collector now accepts only the real matrix
  filenames — `<prefix>.c<digits>.txt` for `-gk 1` (centred, default) and
  `<prefix>.s<digits>.txt` for `-gk 2` (standardized) — so `.log.txt` can never
  be mistaken for a matrix. The stage passes its `-gk` mode so the expected
  suffix is known.
- **`akgwas build` crashed when the kmtricks run directory already existed.**
  The `build` stage created `kmatrixs_out/` before calling `kmtricks pipeline`,
  but kmtricks aborts if `--run-dir` exists (`[CheckFailedError] ... Directory
  already exists!`). The run directory is now cleared and left for kmtricks to
  create itself, which also makes `--force` re-runs (and the smoke test's forced
  second pass) work.
- **Failures were silent unless `--verbose`.** `run_module` caught
  `AkgwasError` and wrote it only to the log file, so a failed run ended with a
  bare non-zero exit code and no message on the console. `Logger.error` now
  always echoes to stderr, and the `kgwas` stage appends the tail of GEMMA's
  per-trait log when `gemma -lmm` exits non-zero, so numerical LMM failures
  ("brent did not converge") are visible instead of hidden behind `E601`.
- **The smoke test's synthetic data made GEMMA's LMM fail.** A 6-sample,
  2-variant dataset yields a rank-deficient GRM: GEMMA printed "Matrix G has 3
  eigenvalues close to zero", the GSL `brent` optimiser failed ("function value
  is not finite"), and the run exited non-zero. (A panel made entirely of *rare*
  k-mers fails the same way -- the GRM is poorly conditioned.) `make_synthetic`
  now plants a wide panel of **common** k-mers (`4 * n_samples` k-mers, each
  carried by 25-75% of the samples) and the smoke test uses 20 samples, so the
  GRM is well conditioned and the LMM converges.
- **The smoke test's synthetic panel still created spurious, trait-correlated
  k-mers.** Inserting each planted k-mer into a random background also creates
  *flanking* k-mers that overlap the insertion. When several carriers happen to
  share the single base adjacent to an insertion, they also share a flanking
  k-mer; that flanking variant is a strict subset of the planted k-mer's
  carriers, so it correlates with the trait and pollutes both the association
  test and the GRM. The failing run showed significant k-mers with `af=0.10`
  (i.e. 2/20 carriers), which the panel design could never produce, while the
  planted causal k-mer fell out of `sigK`. `make_synthetic` now gives each
  **causal** k-mer at most four carriers and writes **distinct** flanking bases
  for each carrier, so each flanking k-mer is unique to one genome and is dropped
  by `recurrence-min` (it occurs only once). Only the causal k-mers need this:
  flankers of the common panel k-mers are trait-independent, so they are
  harmless. Each k-mer is also canonicalised before planting, so the run is
  reproducible regardless of the counter's strand handling.
- The smoke test now accepts the planted k-mer on either strand (a k-mer counter
  may store the canonical form).
- Phenotype/fof/sample order mismatches that silently produced wrong
  associations: the ID set is now checked and rows reordered to `.fam` order.
- First significant k-mer could be dropped by an `NR>1` off-by-one when building
  the FASTA of significant k-mers — the k-mer is now taken directly from the
  sigK row.
- `--multibed-threads` was ignored: the streamed `kmtricks aggregate` step now
  honours it (it previously used `--build-threads`); `--qk-threads` is now
  passed to `plink --pca` (previously ignored); the vestigial `--kgwas-threads`
  flag was dropped (GEMMA jobs are single-threaded — use `--cores`).
- A possible deadlock when `kmtricks aggregate` wrote more than a pipe buffer to
  stderr while its stdout was being consumed — stderr now goes to a temp file.
- `readme` filename case (`Readme.md` → `README.md`) so the build backend finds
  it on case-sensitive filesystems.
- `recipe/meta.yaml` `about.home` pointed at the wrong repository name.

### Removed

- The configuration file entirely (`config.yaml` / any interim `config.toml`):
  all settings are command-line flags now.
- Runtime dependencies `snakemake`, `perl`, `parallel`, `rush`.
- The intermediate multi-gigabyte text matrix and `.tped` files: the matrix is
  streamed from `kmtricks aggregate` and, in the default backend, no `.tped` is
  ever written.

## [0.1.0]

Initial release: a Snakemake workflow driving kmtricks, plink, GEMMA, and
bwa/samtools via Bash scripts. Configuration in `config.yaml`; one
single-column phenotype file per trait under `pheno_dir/`.
