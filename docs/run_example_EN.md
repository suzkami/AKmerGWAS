# AKmerGWAS run examples (v6)

This document provides copy-ready commands: a single `akgwas all` run, the three steps run separately against the same bed, and the external-bed / external K/PCA variants. It ends with the complete parameter list and the required inputs for each step.

> See also: `docs/AKmerGWAS_Q&A.md` (hard requirements for external beds, why v6 cannot simply skip PCA).
> Defaults and behavior descriptions match the v6 source; they may change in later versions.

---

## 1. Prepare the inputs

Example layout:

```text
/data/akgwas_demo/
├── fof.txt                 # sample -> input files
├── pheno.txt               # multi-column phenotype table
├── ref_list.txt            # reference index list
├── reads/                  # reads / assemblies / genomes, any format
│   ├── A88.fq.gz
│   └── A108.fq.gz
└── refs/                   # bowtie or bwa index prefix
    ├── TAIR10.1.bt2
    ├── TAIR10.2.bt2
    └── ...
```

`fof.txt`: one sample per line, `name : file1 ; file2 ; ...`. Lines starting with `#` are comments.

```text
A88  : /data/akgwas_demo/reads/A88.fq.gz
A108 : /data/akgwas_demo/reads/A108.fq.gz
# paired-end example: S1 : /path/S1_R1.fq.gz ; /path/S1_R2.fq.gz
```

`pheno.txt`: first column is the sample ID, the first row is the header (trait names), missing values are `NA`.

```text
ID      ft10    ft16
A88     61.25   22.1
A108    70      NA
```

`ref_list.txt`: one index prefix per line. For bowtie use the `bowtie-build` prefix; for bwa use the `bwa index` prefix (in bwa mode an entry ending in `.sa` is also accepted and stripped automatically).

```text
/data/akgwas_demo/refs/TAIR10
```

The sample identifiers must agree in all three places: names in `fof.txt`, IDs in `pheno.txt`, and IIDs in the bed `.fam`. `--programhome` must be an absolute path; all other relative paths are resolved against it.

---

## 2. Full pipeline in one command: `akgwas all`

`all` runs build → multibed → qk → kgwas → collect → alignment in order.

```bash
akgwas all \
    --programhome /data/akgwas_demo \
    --cores 32 \
    --kmer-len 31 \
    --fof fof.txt \
    --pheno pheno.txt \
    --hard-min 2 \
    --recurrence-min 5 \
    --partitions 0 \
    --build-threads 8 \
    --multibed-threads 8 \
    --block-size 1G \
    --maf 0.01 \
    --bed-backend python \
    --sample-fraction 0.2 \
    --seed 123456789 \
    --gk 1 \
    --qk-threads 4 \
    --threshold 1e-4 \
    --ref-list ref_list.txt \
    --aligner bowtie \
    --alignment-threads 4 \
    --gemma-combine
```

Notes:

- `--cores` controls exactly two pools: the (chunk × trait) jobs in kgwas and the per-reference jobs in alignment. Parallelism inside build / multibed / qk is controlled by `--build-threads`, `--multibed-threads`, and `--qk-threads`.
- To preview without executing anything, add `--dry-run`.
- If a run fails midway, re-running the same command skips completed stages via the markers in `.akgwas/`; add `--force` to re-run them.
- `--gemma-combine` additionally writes the full all-k-mers association dump (`gemma_combine/`, potentially very large). It is off by default.

---

## 3. Three separate steps (same bed)

All three steps must use the same `--programhome` (and the same `--outdir` / `--bed-dir` / `--run-name`, if you changed them). Build-only flags passed to kgwas or alignment are accepted but have no effect.

### Step 1: `akgwas build` (counting + chunking)

```bash
akgwas build \
    --programhome /data/akgwas_demo \
    --cores 32 \
    --kmer-len 31 \
    --fof fof.txt \
    --hard-min 2 \
    --recurrence-min 5 \
    --partitions 0 \
    --build-threads 8 \
    --multibed-threads 8 \
    --block-size 1G \
    --maf 0.01 \
    --bed-backend python
```

This produces `<programhome>/kmatrixs_out/` and `<programhome>/bed/matrix_0001.{bed,bim,fam}`, ..., plus `bed/sample_id_inorder.txt`. These `matrix_*` files are the bed chunks reused by kgwas and alignment.

### Step 2: `akgwas kgwas` (QK + association + collection)

```bash
akgwas kgwas \
    --programhome /data/akgwas_demo \
    --cores 32 \
    --fof fof.txt \
    --pheno pheno.txt \
    --sample-fraction 0.2 \
    --seed 123456789 \
    --gk 1 \
    --qk-threads 4 \
    --threshold 1e-4 \
    --gemma-combine
```

This module contains three internal stages: qk (sample the bed chunks from step 1 to build `QK/QK_kinship.txt` and `QK/QK_pca1-3.txt`), kgwas (one GEMMA call per chunk × trait), and collect (merge sigK, fasta, siggeno, ...).

### Step 3: `akgwas alignment` (map back to references)

```bash
akgwas alignment \
    --programhome /data/akgwas_demo \
    --cores 8 \
    --fof fof.txt \
    --kmer-len 31 \
    --ref-list ref_list.txt \
    --aligner bowtie \
    --alignment-threads 4
```

Important:

- `--kmer-len` must match step 1 (31 by default). If it differs, the CIGAR `<k>M` / `NM:i:0` filter drops every hit and the output is empty.
- In v6, preflight always checks that the fof file exists, so alignment also needs `--fof`. Its contents are not used by alignment, but the file must exist.
- Multiple GWAS runs can share the bed and QK: run kgwas with `--run-name ft10` and `--run-name ft16` to isolate phenotype-dependent outputs under `<outdir>/<run-name>/`; alignment must use the same `--run-name`.

---

## 4. External bed / external K and PCA

### 4.1 External K/PCA only (bed still built by akgwas)

`akgwas all` works directly with external K/PCA:

```bash
akgwas all \
    --programhome /data/akgwas_demo \
    --cores 32 \
    --fof fof.txt \
    --pheno pheno.txt \
    --recurrence-min 5 \
    --maf 0.01 \
    --use-external-qk \
    --kinship K.txt \
    --pca PCA.txt \
    --pca-has-id \
    --threshold 1e-4 \
    --ref-list ref_list.txt \
    --aligner bowtie
```

When running step by step, add the external arguments only to the kgwas step:

```bash
akgwas build --programhome /data/akgwas_demo --cores 32 \
    --fof fof.txt --recurrence-min 5 --maf 0.01

akgwas kgwas --programhome /data/akgwas_demo --cores 32 \
    --fof fof.txt --pheno pheno.txt \
    --use-external-qk --kinship K.txt --pca PCA.txt --pca-has-id

akgwas alignment --programhome /data/akgwas_demo --cores 8 \
    --fof fof.txt --kmer-len 31 --ref-list ref_list.txt
```

Formats:

- **K (kinship)**: GEMMA `-k` format, first line is the sample count n, followed by n rows of n values. Row/column order must match the bed fam order; v6 does not verify sample IDs in K.
- **PCA / covariates**: with `--pca-has-id`, sample IDs are taken from column 2 and columns 3–5 are used as covariates (at most 3). Without it, the file is assumed headerless with numeric columns, row order equals the bed fam order, and only the row count is checked.

### 4.2 External bed + external K/PCA (do not run `all`)

Once the external bed is in `<outdir>/bed/`, **do not run `build`, `multibed`, or `all`**: multibed deletes and rebuilds `bed/matrix_*`. The downstream equivalent of `all` is the kgwas module plus alignment:

```bash
akgwas kgwas \
    --programhome /data/akgwas_demo \
    --cores 32 \
    --fof fof.txt \
    --pheno pheno.txt \
    --use-external-qk \
    --kinship K.txt \
    --pca PCA.txt \
    --pca-has-id \
    --threshold 1e-4

akgwas alignment \
    --programhome /data/akgwas_demo \
    --cores 8 \
    --fof fof.txt \
    --kmer-len 31 \
    --ref-list ref_list.txt \
    --aligner bowtie
```

To chain them into one command (`&&` ensures alignment only runs after kgwas succeeds):

```bash
akgwas kgwas --programhome /data/akgwas_demo --cores 32 \
    --fof fof.txt --pheno pheno.txt \
    --use-external-qk --kinship K.txt --pca PCA.txt --pca-has-id \
  && akgwas alignment --programhome /data/akgwas_demo --cores 8 \
    --fof fof.txt --kmer-len 31 --ref-list ref_list.txt
```

Hard requirements for an external bed (details in `docs/AKmerGWAS_Q&A.md`):

1. Files are `<outdir>/bed/matrix_XXXX.bed/.bim/.fam`, PLINK1 SNP-major (any `plink --make-bed` output works).
2. All chunks share exactly the same `.fam` sample order; writing `bed/sample_id_inorder.txt` (one IID per line, in fam order) is recommended.
3. **Column 2 of `.bim` must be the k-mer sequence itself**; otherwise collect finds no variants and `siggeno` / alignment output is empty.
4. `fof.txt` must exist and its sample names must match `pheno.txt` and the fam IIDs (when build is not run, the file paths inside fof may be placeholders).
5. If you reuse an old outdir, add `--force` or use a fresh outdir, so completion markers do not skip qk / kgwas.

### 4.3 External bed + internal QK

Simply omit `--use-external-qk`; qk samples the external bed and rebuilds K/PCA:

```bash
akgwas kgwas \
    --programhome /data/akgwas_demo \
    --cores 32 \
    --fof fof.txt \
    --pheno pheno.txt \
    --sample-fraction 0.2 \
    --seed 123456789 \
    --gk 1 \
    --qk-threads 4

akgwas alignment \
    --programhome /data/akgwas_demo \
    --cores 8 \
    --fof fof.txt \
    --kmer-len 31 \
    --ref-list ref_list.txt
```

---

## 5. Complete parameter list

Note: all four subcommands (`build`, `kgwas`, `alignment`, `all`) share the same argument definitions. **At the argument level, `--programhome` is the only required one**; everything else falls back to its default. "Applies to" shows where a flag actually has an effect; passing it to another step is accepted but ignored.

### 5.1 Common / input-output

| Flag | Default | Applies to | Required? | Description |
| --- | --- | --- | --- | --- |
| `--programhome` | none | all | **required** | Absolute base directory for paths and outputs |
| `--cores` | `1` | kgwas, alignment | optional | Parallel jobs: (chunk × trait) for kgwas, references for alignment |
| `--dry-run` | off | all | optional | Print commands without writing anything |
| `--force` | off | all | optional | Ignore completion markers and re-run |
| `-q`, `--quiet` | off | all | optional | Less console output (errors still printed) |
| `--outdir` | programhome | all | optional | Output root; relative paths are resolved against programhome |
| `--run-name` | empty | kgwas, collect, alignment | optional | Isolate phenotype-dependent outputs under `<outdir>/<run-name>/`; bed and QK stay shared |
| `--fof` | `fof.txt` | build, kgwas; alignment only checks existence | optional | Sample list file |
| `--pheno` | `pheno.txt` | kgwas | optional | Multi-column phenotype table |
| `--kmer-len` | `31` | build, alignment | optional | k-mer length; alignment must match the build |
| `--ulimit-min` | `4096` | all (preflight) | optional | Minimum open-file limit |
| `--bed-dir` | `bed` | build, qk, kgwas, collect | optional | Shared bed chunk directory name |
| `--qk-dir` | `QK` | qk, kgwas | optional | Shared K/PCA directory name |
| `--temp-dir` | `temp_gemma` | kgwas, collect | optional | Per-run GEMMA work directory name |
| `--final-dir` | `kmergwas_outputs` | collect, alignment | optional | Per-run results directory name |
| `--kmtricks-run-dir` | `kmatrixs_out` | build, multibed | optional | kmtricks run directory name |

### 5.2 build (k-mer counting)

| Flag | Default | Applies to | Required? | Description |
| --- | --- | --- | --- | --- |
| `--hard-min` | `2` | build | optional | Minimum within-sample count to call a k-mer present |
| `--recurrence-min` | `2` | build | optional | Minimum cross-sample recurrence; passed to both kmtricks count and merge |
| `--partitions` | `0` | build | optional | kmtricks partition count; `0` = automatic |
| `--build-threads` | `8` | build | optional | kmtricks count threads (merge is fixed at 2) |
| `--aggregate-to-stdout`, `--no-aggregate-to-stdout` | on | build / multibed | optional | multibed consumes aggregate output from stdout; when off, a temporary file is written first |

### 5.3 multibed (chunking, part of the build module)

| Flag | Default | Applies to | Required? | Description |
| --- | --- | --- | --- | --- |
| `--block-size` | `1G` | multibed | optional | Target byte budget per chunk; binary units (K/M/G/T) |
| `--min-block-size` | `16M` | multibed | optional | Only used for validation (requires min ≤ block-size ≤ max); it does not clamp or merge chunks at runtime |
| `--max-block-size` | `10G` | multibed | optional | Same as above |
| `--maf` | `0.01` | multibed | optional | MAF filter: keep if `max(#absent, #present) <= int(n - maf*n)` |
| `--multibed-threads` | `8` | multibed | optional | `kmtricks aggregate` threads |
| `--bed-backend` | `python` | multibed | optional | `python` writes bed directly; `plink` writes tped first and converts with plink |

### 5.4 qk (kinship and PCA, part of the kgwas module)

| Flag | Default | Applies to | Required? | Description |
| --- | --- | --- | --- | --- |
| `--use-external-qk` | off | qk | optional | Use external K/PCA; requires both `--kinship` and `--pca` |
| `--kinship` | empty | qk | conditionally required | External kinship matrix (relative to programhome) |
| `--pca` | empty | qk | conditionally required | External PCA / covariate file (relative to programhome) |
| `--sample-fraction` | `0.2` | qk | optional | Fraction of chunks sampled for internal K/PCA |
| `--seed` | `123456789` | qk | optional | Sampling seed (reproducibility) |
| `--gk` | `1` | qk | optional | GEMMA `-gk`: `1` centred, `2` standardized |
| `--pca-has-id` | off | qk | optional | External PCA carries sample IDs (column 2); validate and reorder by ID |
| `--qk-threads` | `1` | qk | optional | `plink --pca` threads |

### 5.5 kgwas (association)

| Flag | Default | Applies to | Required? | Description |
| --- | --- | --- | --- | --- |
| `--threshold` | `0.0001` | kgwas | optional | Significance cutoff: k-mers with `p_wald <= threshold` go to sigK |

### 5.6 collect (merging)

| Flag | Default | Applies to | Required? | Description |
| --- | --- | --- | --- | --- |
| `--gemma-combine`, `--no-gemma-combine` | off | collect | optional | Also write the full all-k-mers association dump (`gemma_combine/`, potentially very large) |

### 5.7 alignment (mapping to references)

| Flag | Default | Applies to | Required? | Description |
| --- | --- | --- | --- | --- |
| `--ref-list` | `ref_list.txt` | alignment | optional | Reference index list (relative to programhome) |
| `--aligner` | `bowtie` | alignment | optional | `bowtie` or `bwa` |
| `--max-hits` | `5000` | alignment | optional | Drop an entire k-mer when it has more hits than this |
| `--alignment-threads` | `1` | alignment | optional | Threads per aligner process (bowtie `-p` / bwa `-t`) |

### 5.8 Top level

| Flag | Description |
| --- | --- |
| `akgwas --version` | Print the version |
| `akgwas --errors` | Print the error-code table |
| `akgwas --help` / `akgwas <subcommand> --help` | Help |

---

## 6. Required inputs per step (quick reference)

Strictly speaking, the only required argument is `--programhome`. The "required data" column lists what must exist before that step can run; files may sit at their default paths or be pointed to with the corresponding flag.

| Step | Required arguments | Required data / tools | Outputs |
| --- | --- | --- | --- |
| `build` | `--programhome` | The file given by `--fof` (and every input file it references); `kmtricks`; `plink` only when `--bed-backend plink` | `kmatrixs_out/`, `bed/matrix_*`, `bed/sample_id_inorder.txt` |
| `kgwas` | `--programhome` | `--fof` and `--pheno` files; `bed/matrix_*.bed` (from build or external); `gemma` is always needed, internal QK additionally needs `plink`, external QK additionally needs the `--kinship` and `--pca` files | `QK/`, `<run_root>/temp_gemma/`, `<run_root>/kmergwas_outputs/` |
| `alignment` | `--programhome` | The `--fof` file (preflight only); `--ref-list`; `kmergwas_outputs/sigk/*_sigK.txt` and `fasta/` from the previous step; `bowtie` or `bwa` | `kmergwas_outputs/samfile/`, `kmergwas_outputs/sigk_mapped/` |
| `all` | `--programhome` | All of the above (`fof`, `pheno`, `ref_list`, plus kmtricks, gemma, plink, and the aligner) | All of the above |

---

## 7. Re-runs, skipping, and `--force`

- After each stage succeeds, a marker is written under `<outdir>/.akgwas/`; running the same command again skips completed stages.
- Markers are based on input timestamps/sizes and **do not include a parameter signature**. After changing result-affecting parameters (`--hard-min`, `--recurrence-min`, `--maf`, `--kmer-len`, ...), add `--force`; if necessary, delete `kmatrixs_out/`, `bed/`, and `QK/` first so old and new results cannot mix.
- To re-run just one step, call that subcommand with `--force`.

---

## 8. Where the outputs go

| Path | Content |
| --- | --- |
| `<outdir>/kmatrixs_out/` | kmtricks counting results |
| `<outdir>/bed/matrix_XXXX.{bed,bim,fam}` | Bed chunks used by kgwas / alignment |
| `<outdir>/bed/sample_id_inorder.txt` | Sample order |
| `<outdir>/QK/QK_kinship.txt` | Kinship matrix (stable name) |
| `<outdir>/QK/QK_pca1-3.txt` | First 3 PCs (stable name) |
| `<outdir>/QK/QK_samples.txt` | Sample order used by qk |
| `<run_root>/pheno_single/` | One single-column phenotype file per trait |
| `<run_root>/temp_gemma/` | Per-chunk GEMMA intermediate results |
| `<run_root>/kmergwas_outputs/sigk/` | Significant k-mers per trait (`<trait>_sigK.txt`) |
| `<run_root>/kmergwas_outputs/fasta/` | FASTA of significant k-mers for alignment |
| `<run_root>/kmergwas_outputs/siggeno/` | Genotype bed of significant k-mers |
| `<run_root>/kmergwas_outputs/samfile/` | Alignment SAM files |
| `<run_root>/kmergwas_outputs/sigk_mapped/` | Results with reference positions (`<trait>_<ref>_chrpos_sigk.txt`) |
| `<outdir>/.akgwas/` | Completion markers, logs, run manifest |

`<run_root>` is `<outdir>` unless `--run-name` is given, in which case it is `<outdir>/<run-name>`.

---

## 9. Troubleshooting

- **`E601` / `QK_kinship.cXX.txt` not found**: in GEMMA's real output filename, `XX` is literal. v6 accepts `cXX` / `sXX`; if it still fails, check the actual names under `QK/output/` and your GEMMA version.
- **`no bed chunks under .../bed`**: `--outdir` / `--bed-dir` differs from the previous step, or the external bed is not named `matrix_*.bed`.
- **Sample ID mismatch**: fof names, pheno IDs, and bed fam IIDs must agree; for an external bed, check the IID column of the fam first.
- **alignment output is empty**: `--kmer-len` differs from the build, or there are no significant k-mers (check whether `sigk_mapped/` has files and whether `samfile/` is empty).
- **External bed produces no siggeno**: column 2 of `.bim` is not the k-mer sequence (see `docs/AKmerGWAS_Q&A.md`).
