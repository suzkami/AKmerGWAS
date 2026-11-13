# AKmerGWAS

A k-mer based genome-wide association study (GWAS) pipeline that accepts any
input format — short reads, assemblies, or whole genomes.

The pipeline counts k-mers with [kmtricks](https://github.com/tlemane/kmtricks),
turns the presence/absence matrix into PLINK bed chunks, runs
[GEMMA](https://github.com/genetics-statistics/GEMMA) association per chunk and
trait, and maps significant k-mers back to reference genomes with
`bowtie` or `bwa`.

Version **1.0.0** is a ground-up rewrite of the Python/Snakemake/Bash stack into
a single pure-Python orchestrator (standard library only), driven **entirely
from the command line** — there is no config file. See
[Migration from 0.1.0](#migration-from-010).

---

## Installation

```bash
conda install -c bioconda akgwas
akgwas --version
```

Runtime external tools (declared as conda dependencies):
`kmtricks`, `gemma`, `plink`, `bowtie`, `bwa`.

---

## Quick start

The pipeline maps the significant k-mers back to reference genomes, so it needs
a **reference index**. akgwas never builds the index itself — with the default
aligner (`bowtie`), build one once per reference with `bowtie-build` and write
the resulting index **prefix** into `ref_list.txt`:

```bash
# 1. build a bowtie index (one-off; reuse across runs)
mkdir -p /data/indexes
bowtie-build /data/genomes/pseudo10001.fasta /data/indexes/pseudo10001
# -> writes /data/indexes/pseudo10001.{1,2,3,4}.ebwt and .rev.{1,2}.ebwt

# 2. ref_list.txt: one index prefix per line -- the path you passed to
#    bowtie-build, WITHOUT the .1.ebwt / .rev.1.ebwt suffix
echo /data/indexes/pseudo10001 > /data/run/ref_list.txt
```

(With `--aligner bwa`, list bwa index prefixes instead — see the
`ref_list.txt` section under [Inputs](#inputs).)

Then put `fof.txt`, `pheno.txt` and `ref_list.txt` in `--programhome` and run the
whole pipeline:

```bash
# Required files in --programhome: fof.txt, pheno.txt, ref_list.txt
akgwas all --programhome /data/run --cores 32
```

That is the whole interface. Every value has a default; pass the matching flag
to override it. Only `--programhome` (the base for relative paths and all
outputs) is required.

Run individual phases:

```bash
akgwas build     --programhome /data/run --cores 8    # kmtricks count/merge + bed chunking
akgwas kgwas     --programhome /data/run --cores 32   # K/PCA + GEMMA + collect
akgwas alignment --programhome /data/run --cores 8    # map significant k-mers
```

Preview without running anything, or re-run ignoring saved progress:

```bash
akgwas all --programhome /data/run --dry-run
akgwas all --programhome /data/run --force
```

Run a subcommand with `--help` to see the grouped options (input/output, build,
multibed, qk, kgwas, alignment, collect), e.g.
`akgwas build --programhome /data/run --block-size 5G --maf 0.05`.

---

## Inputs

Three plain-text files live in `--programhome` (names overridable with `--fof`,
`--pheno`, `--ref-list`).

### 1. `fof.txt` — sequences

One sample per line: `name : file1 ; file2 ; ...`.

```
A88  : /data/pseudo88_1.fasta.gz ; /data/pseudo88_2.fasta.gz
A108 : /data/pseudo108.fasta.gz
```

### 2. `pheno.txt` — phenotype table

A single table. The **first column holds sample IDs**, the **header row holds
trait names**, and missing values are written as `NA`.

```
ID      ft10     ft16
A88     61.25    22.1
A108    70       NA
A139    68.75    19.7
```

For GWAS the pipeline: verifies that the sample-ID *set* equals the bed `.fam`
set, reorders rows to the `.fam` order, and exports one single-column,
header-less file per trait (`<run_root>/pheno_single/<trait>.txt`) for GEMMA `-p`.

### 3. `ref_list.txt` — reference indices

One entry per line. With `--aligner bowtie` (default) entries are **bowtie index
prefixes** — the `<prefix>` you pass to `bowtie-build <ref.fasta> <prefix>`.
akgwas does **not** run `bowtie-build`; prepare the index first, then write its
path here:

```bash
bowtie-build /data/genomes/pseudo10001.fasta /data/indexes/pseudo10001
```

```
# ref_list.txt
/data/indexes/pseudo10001
/data/indexes/pseudo10002
```

With `--aligner bwa` entries are bwa index prefixes (or `*.sa` paths).

---

## The command line is the config

Every configuration value is a flag with a default; nothing is read from a file.
The exact invocation is recorded in the run manifest, so a run stays
reproducible. Highlights:

| Flag | Default | Meaning |
|---|---|---|
| `--programhome` | *(required)* | base dir for paths and outputs |
| `--outdir` | `programhome` | output root (relative → programhome) |
| `--run-name` | *(none)* | isolate per-phenotype outputs, see below |
| `--kmer-len` | `31` | k-mer length |
| `--hard-min` | `2` | min within-sample count to call present (use `1` for assemblies/genomes) |
| `--block-size` / `--min-block-size` / `--max-block-size` | `1G` / `16M` / `10G` | chunk size bounds |
| `--maf` | `0.01` | minor-allele-frequency filter |
| `--bed-backend` | `python` | `python` (direct bed) or `plink` |
| `--use-external-qk` `--kinship` `--pca` | off | supply external K/PCA |
| `--sample-fraction` `--seed` `--gk` | `0.2` / `123456789` / `1` | internal K/PCA |
| `--threshold` | `0.0001` | `p_wald` cutoff for extraction |
| `--aligner` | `bowtie` | `bowtie` or `bwa` |
| `--max-hits` | `5000` | drop a k-mer with more hits than this |
| `--gemma-combine` | off | also write the full all-k-mers dump (`gemma_combine/`; can be very large) |

Run `akgwas <subcommand> --help` for the complete list.

### Reusing bed/QK across phenotype runs — `--run-name`

`build`, `multibed` and `qk` are **trait-independent**, so their outputs
(`bed/`, `QK/`, `kmatrixs_out/`) live directly under `outdir` and are shared.
Everything that depends on the phenotype (`pheno_single/`, `temp_gemma/`,
`kmergwas_outputs/`, and the kgwas/collect/alignment completion markers) lives
under `run_root` = `outdir` (no run name) or `outdir/<run-name>`.

That lets you compute the expensive matrix + K/PCA once and run several
phenotype sets at different times without clobbering each other or wrongly
skipping stages:

```bash
# once: matrix + K/PCA + first trait set
akgwas all  --programhome /data/run --cores 32

# later: a different phenotype table, reusing bed/ and QK/
akgwas kgwas     --programhome /data/run --run-name ft16 --pheno pheno_ft16.txt
akgwas alignment --programhome /data/run --run-name ft16
```

`ft16` lands under `/data/run/ft16/kmergwas_outputs/`, with its own marker file
so it does not collide with (or get skipped because of) the first run.

---

## Sample alignment — the one invariant that matters

Everything is conformed to **the bed `.fam` order**, which is itself derived from
`fof.txt` (sample declaration order — *not* sorted):

* the phenotype table is reordered to the `.fam` order (ID-based match);
* the kinship/PCA matrices must be built with **the same sample order as the
  final GWAS bed chunk**;
* the tool enforces this: the internal k-mer K/PCA is built by concatenating bed
  chunks (never `plink --merge-list`, which reorders samples lexicographically),
  and the PCA `eigenvec` IDs are checked against the canonical order.

**When you supply external K/PCA** (`--use-external-qk`), you are responsible for
the sample order. If your PCA file carries sample IDs, add `--pca-has-id` and the
tool will reorder/verify it; the kinship matrix (which has no IDs) is assumed to
follow the same order — see `docs/troubleshooting.md`.

---

## Outputs

```
<outdir>/
├── bed/                     matrix_0001.{bed,bim,fam} ... (trait-independent)
├── QK/                      QK_kinship.txt, QK_pca1-3.txt, QK_samples.txt
├── .akgwas/                 akgwas.log, run_manifest[_<run>].json, <stage>.done
└── <run_root>/              <run_root> = . (no --run-name) or <run-name>/
    ├── pheno_single/        <trait>.txt (single column, fed to GEMMA)
    └── kmergwas_outputs/
        ├── sigk/            <trait>_sigK.txt          (k-mer af beta se p_wald)
        ├── siggeno/         <trait>_bed.{bed,bim,fam} (PLINK)
        ├── fasta/           <trait>_fasta.txt         (sequences for alignment)
        ├── gemma_combine/   <trait>.assoc_combined.txt  (only with --gemma-combine)
        ├── logs/            <trait>.log_combined.txt
        ├── samfile/         <trait>_<ref>.sam
        └── sigk_mapped/     <trait>_<ref>_chrpos_sigk.txt
```

The `build` phase keeps only the kmtricks run directory and the bed chunks: the
multi-gigabyte text matrix and all `.tped` intermediates are never written to
disk (streamed / one chunk at a time).

---

## Resources, disk and reproducibility

* **Disk.** The dominant cost is the kmtricks run directory plus the final bed.
  Chunk size (`--block-size`) bounds peak transient disk.
* **RAM / CPU.** k-mer counting is bounded by the chunk size; GEMMA runs one
  thread per job and concurrency comes from `--cores`.
* **Reproducibility.** The chunk-sampling seed is `--seed`; the tool version,
  the full resolved configuration and selected tool versions are written to
  `.akgwas/run_manifest[_<run-name>].json`.
* **Resume.** Each stage writes an atomic `<stage>.done` marker (with an input
  fingerprint). Re-running skips stages whose inputs are unchanged; `--force`
  ignores markers. Trait-dependent stages are marked per `--run-name`.

---

## Migration from 0.1.0

| Area | 0.1.0 | 1.0.0 |
|---|---|---|
| Config | `config.yaml` | **none — command-line flags** |
| Phenotype | one single-column file per trait under `pheno_dir/` | one multi-column `pheno.txt` (ID + trait columns) |
| Orchestration | Snakemake (`.smk`) + GNU parallel/rush + Bash | pure Python |
| Aligner default | `bwa mem` | `bowtie -a -v 0` (bwa `-k31 -T31 -a -c 5000` available) |
| CLI | `--build/--kgwas/...` flags | `akgwas build\|kgwas\|alignment\|all` |
| Shared vs per-run | all outputs in one tree | `bed/`+`QK/` shared; per-trait outputs under `--run-name` |
| Dropped deps | snakemake, perl, rush, parallel | — |

---

## Development / tests

```bash
pip install -e .
pytest tests/              # unit tests run without external tools
pytest tests/test_smoke.py # end-to-end, needs kmtricks/gemma/plink/bowtie
python -m tests.make_synthetic /tmp/syn   # inspect a generated dataset
```
