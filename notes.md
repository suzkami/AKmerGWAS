# AKmerGWAS v6 — Usage & Troubleshooting Notes

This note records the key operational points for v6: how to recover when `build`
fails, a few usage clarifications (which stage builds QK / `--run-name`
isolation / running bowtie and bwa separately / `--kmer-len`), and one
command-line limitation (subcommands cannot be stacked).

---

## 1. How to re-run when build fails

If `build` fails, or if you are unsure the artifacts are in a good state,
**delete the `kmatrixs_out/`, `bed/` and `QK/` directories, then re-run with
`--force`**:

```bash
rm -rf /data/run/kmatrixs_out /data/run/bed /data/run/QK
akgwas build --programhome /data/run --cores 8 --force
```

Notes:

- Before calling `kmtricks pipeline`, the `build` stage itself does
  `remove_if_exists(kmatrixs_out/)` (kmtricks refuses an existing `--run-dir`),
  so `kmatrixs_out/` usually clears itself; the directories you actually need to
  delete by hand are the **leftover `bed/` and `QK/`**.
- `--force` ignores completion markers and re-runs every stage.
- Why do "delete artifacts + `--force`" together: a completion marker only
  fingerprints the **inputs** — it does **not** check whether the artifacts still
  exist, nor whether parameters changed. So if `fof.txt` is unchanged but the
  artifacts were deleted, the stage misreads it as `valid -> skip` and then dies
  on the missing artifact (the classic pattern is `build: marker valid -> skip`
  immediately followed by `[E401] kmtricks aggregate ... [--run-dir
  .../kmatrixs_out] Directory doesn't exist!`). Deleting the artifacts and adding
  `--force` resets that state cleanly.

---

## 2. Usage clarifications

### 2.1 Which stage builds QK? Does running build alone build QK?

Stage/module map (`MODULE_STAGES` in `pipeline.py`):

| Module | Stages |
|---|---|
| `build` | `build` + `multibed` |
| `kgwas` | `qk` + `kgwas` + `collect` |
| `alignment` | `alignment` |
| `all` | all six |

- `akgwas build` = kmtricks count + merge (produces `kmatrixs_out/`) **then**
  multibed (`kmtricks aggregate` streamed into `bed/matrix_XXXX.{bed,bim,fam}`).
  **It does not build QK.**
- **QK is created in the `qk` stage**, and `qk` is the **first** stage of `kgwas`
  (and of `all`). Its outputs: `QK/QK_kinship.txt`, `QK/QK_pca1-3.txt`,
  `QK/QK_samples.txt` (the intermediate `QK_bed.*` is deleted once used).
- QK is trait-independent and **shared** (under `outdir/QK`), and its completion
  marker is not run-name-scoped, so it is **built only once** and reused across
  `--run-name`s.

### 2.2 How does each phenotype start a new run? Where is the run dir? What about temp_gemma / kmergwas_outputs?

Isolate with `--run-name`:

```bash
akgwas kgwas     --programhome /data/run --run-name ft16 --pheno pheno_ft16.txt --cores 32
akgwas alignment --programhome /data/run --run-name ft16 --cores 8
```

- `--pheno` resolves relative to `programhome` (or give an absolute path).
- `run_root = outdir/<run-name>` (`outdir` defaults to `programhome`), so this
  run's directory is **`/data/run/ft16/`**, containing:
  - `pheno_single/` — one single-column file per trait, fed to GEMMA
  - `temp_gemma/` — one subdir per chunk: `output/<trait>.assoc.txt`, `sigk/`, `fasta/`
  - `kmergwas_outputs/` — the final collect + alignment outputs
- The shared dirs `kmatrixs_out/`, `bed/`, `QK/` stay at the top level of
  `/data/run/` and are not recomputed.
- Running a second phenotype:
  - a **different** `--run-name` → `temp_gemma` / `kmergwas_outputs` live under
    separate run_roots, fully isolated (recommended);
  - **reusing the same run-name** and only changing `--pheno` → both share the
    same `temp_gemma` / `kmergwas_outputs`; kgwas re-runs (the pheno path enters
    the input fingerprint, so the fingerprint changes), files are overwritten by
    name (atomic replace), but **trait files left over from the old phenotype
    that the new one does not have are not cleaned up**. So give each phenotype
    its own `--run-name`.

### 2.3 Running bowtie and bwa separately for the same phenotype; do I only change ref_list.txt? Is ksize fixed?

Alignment outputs all live under `<run_root>/kmergwas_outputs/`:

```
samfile/<trait>_<refname>.sam
sigk_mapped/<trait>_<refname>_chrpos_sigk.txt
```

`<refname>` = the basename of the index path in `ref_list.txt`, with the sequence
extension stripped (`.fasta/.fa/.fasta.gz/.gz/.sa` …); duplicate basenames get
`_1`, `_2`.

- **There is no aligner tag in the filename.** If the bowtie and bwa indices share
  a basename (e.g. both called `pseudo10001`), the two alignments'
  `samfile`/`sigk_mapped` **overwrite each other**. To keep both sets, use
  **different `--run-name`s** (e.g. `ft16_bowtie` / `ft16_bwa`), or name the two
  index prefixes differently.
- Switching to bwa is **not just a `ref_list.txt` edit**: also pass
  `--aligner bwa`, and change the *contents* of `ref_list.txt` to **bwa indices**
  (bwa index prefixes, or `*.sa` paths).
- The `alignment` completion marker's inputs **do not include the aligner** — only
  the `ref_list` file fingerprint. If you change `--aligner` but `ref_list.txt`'s
  contents are unchanged, the stage may `valid -> skip`; pass `--force` in that
  case. Normally you switch index paths, so the contents change and it re-runs
  automatically.
- **ksize is neither fixed nor auto-detected from the matrix** — it is the CLI
  `--kmer-len` (default 31). `build` uses it as `kmtricks --kmer-size`;
  `alignment` uses it twice: for the exact-match CIGAR length `<k>M`, and for
  bwa's `-k <k> -T <k>`. So **if build used a non-default value, you must pass
  the same `--kmer-len` when running alignment alone**; otherwise (e.g. build
  used 25, alignment defaults to 31) every hit is filtered out, `sigk_mapped` is
  empty, and bwa's `-k/-T` are wrong too. (`kgwas`/`collect`/`qk`/`multibed` do
  not use `kmer_len`.)

---

## 3. Subcommands cannot be stacked

A command takes exactly **one** subcommand (argparse subparsers).

- `akgwas kgwas alignment ...` → `unrecognized arguments: alignment`.
- `akgwas --programhome /data/run ...` (no subcommand) → it just prints help and
  does nothing.
- To run kgwas then alignment in one go, three options:
  1. two separate commands;
  2. shell chaining:
     ```bash
     akgwas kgwas --programhome /data/run --run-name ft16 --pheno pheno_ft16.txt --cores 32 \
       && akgwas alignment --programhome /data/run --run-name ft16 --cores 8
     ```
  3. use `all` (runs `build→multibed→qk→kgwas→collect→alignment`; the already-done
     `build`/`multibed`/`qk` are skipped via completion markers, so the extra cost
     is just kgwas+collect+alignment):
     ```bash
     akgwas all --programhome /data/run --run-name ft16 --pheno pheno_ft16.txt --cores 32
     ```
- To give kgwas and alignment **different** `--cores`, use two separate commands
  (`all` only takes one `--cores`).
