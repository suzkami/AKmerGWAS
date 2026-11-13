# Troubleshooting

Every failure the pipeline can attribute prints a stable code, e.g.

```
[E301] phenotype sample IDs do not match the bed .fam sample set
    hint: Phenotype/fof sample IDs must form the same set as the bed .fam IDs.
```

Print the whole table at any time with `akgwas --errors`. Codes are grouped by
stage:

| Range | Area |
|---|---|
| `E1xx` | environment / preflight (tools, ulimit, paths) |
| `E2xx` | input validation (fof, phenotype, reference list, config) |
| `E3xx` | sample alignment / sample-order invariants |
| `E4xx` | build + multibed (kmtricks, chunking, bed writing) |
| `E5xx` | kinship / PCA (plink, gemma -gk) |
| `E6xx` | GWAS (gemma -lmm) |
| `E7xx` | collection / merge |
| `E8xx` | alignment (bowtie / bwa / SAM parsing) |
| `E9xx` | internal / unexpected |

---

## Error-code reference

### E1xx — environment / preflight

| Code | Meaning | Fix |
|---|---|---|
| `E101` | a required external tool is not on `PATH` | activate the conda env, or `conda install -c bioconda <tool>` |
| `E102` | open-file limit (`ulimit -n`) below `--ulimit-min` | `ulimit -n 4096` before running, or lower `--ulimit-min` |
| `E103` | working/output directory is not writable | check permissions and free disk space |

### E2xx — input validation

| Code | Meaning | Fix |
|---|---|---|
| `E201` | an input file does not exist | check the flag/path; relative paths resolve against `--programhome` |
| `E202` | invalid setting (bad value, bad block-size range, unknown choice) | check the value; `akgwas <subcommand> --help` lists valid ranges |
| `E204` | `fof.txt` is malformed | each line must be `name : file1 ; file2 ; ...` |
| `E205` | row count mismatch (phenotype vs fof) | one phenotype row per fof sample |
| `E206` | reference list is empty | add at least one index line to `ref_list.txt` |
| `E207` | phenotype table malformed (ragged rows, no header) | first row = header, first column = sample ID |

### E3xx — sample alignment (the invariant)

| Code | Meaning | Fix |
|---|---|---|
| `E301` | phenotype/fof sample-ID set ≠ bed `.fam` ID set | align the sample sets; check for typos, whitespace, duplicate IDs |
| `E302` | K/PCA sample order ≠ GWAS bed order | rebuild K/PCA on the same sample order (see below) |
| `E303` | sample order could not be verified (PCA has no IDs) | pass `--pca-has-id` if the PCA file carries IDs |

### E4xx — build + multibed

| Code | Meaning | Fix |
|---|---|---|
| `E401` | `kmtricks` count/merge failed | read the captured kmtricks log in `.akgwas/` |
| `E402` | kmtricks produced an empty matrix | inputs empty, or `--hard-min`/`--recurrence-min` too strict |
| `E403` | unknown bed backend | use `--bed-backend python` or `--bed-backend plink` |
| `E404` | a matrix row was malformed (wrong width / non-0-1) | usually corrupt intermediate input; re-run `build` with `--force` |

### E5xx — kinship / PCA

| Code | Meaning | Fix |
|---|---|---|
| `E501` | `plink --pca` failed | inspect the plink log; often a malformed bed |
| `E502` | `gemma -gk` failed | inspect the gemma log; usually a kinship input issue |

### E6xx — GWAS

| Code | Meaning | Fix |
|---|---|---|
| `E601` | `gemma -lmm` failed for a chunk×trait | inspect the per-trait gemma log under `logs/` (the tail is also appended to the error) |
| `E602` | a GEMMA `.assoc.txt` had an unexpected column count | check the GEMMA version; expected 12 columns |

> **Numerical LMM failures** ("Matrix G has N eigenvalues close to zero", "GSL
> ERROR: ... brent.c", "Brent did not converge") are reported as `E601`. They
> mean the genomic relationship matrix is rank-deficient — almost always too few
> samples and/or too few k-mer variants (e.g. a tiny pilot dataset). Widen the
> panel (more samples, lower `--hard-min`/`--recurrence-min`, lower `--maf`) or
> supply an external well-conditioned kinship via `--use-external-qk`.

### E7xx — collection

| Code | Meaning | Fix |
|---|---|---|
| `E701` | nothing to collect (no chunks, or no significant k-mers) | for zero hits this is informational; lower `--threshold` |

### E8xx — alignment

| Code | Meaning | Fix |
|---|---|---|
| `E801` | unknown aligner | use `--aligner bowtie` or `--aligner bwa` |
| `E802` | the aligner failed | inspect the aligner stderr; check the index prefix in `ref_list.txt` |

### E9xx — internal

| Code | Meaning | Fix |
|---|---|---|
| `E999` | unexpected internal error | report it together with `.akgwas/run_manifest.json` |

---

## The sample-order invariant (read this if you see E301/E302/E303)

Everything is conformed to **the bed `.fam` order**, and that order comes from
`fof.txt` (sample declaration order — *not* sorted). Three inputs must agree:

1. **Phenotype** — matched to `.fam` by ID, then reordered automatically. If the
   ID *sets* differ you get `E301`.
2. **Kinship / PCA** — must be built on the *same* sample order as the final
   GWAS bed. Without `--use-external-qk`, the tool builds this itself by
   concatenating bed chunks (never `plink --merge-list`, which silently reorders
   samples lexicographically) and verifies the PCA IDs, so you cannot get this
   wrong internally.
3. **External K/PCA** — with `--use-external-qk`, *you* own the order.
   - If your PCA file has sample IDs in its first columns, add
     `--pca-has-id`; the tool reorders and verifies it.
   - The kinship matrix carries **no IDs**, so it is assumed to follow the same
     order you gave the PCA. If in doubt, regenerate K and PCA together from the
     same SNP set with the same sample order.

The GWAS bed `.fam` is the source of truth. To print it:

```bash
cut -d' ' -f2 bed/matrix_0001.fam
```

Compare that column against your external PCA IDs — they must match exactly.

---

## Reproducibility and resume

- Each stage writes an atomic `<stage>.done` marker (with an input fingerprint)
  under `.akgwas/`. A re-run skips stages whose inputs are unchanged.
- `--force` ignores the markers and re-runs everything.
- Trait-dependent stages (`kgwas`, `collect`, `alignment`) are marked per
  `--run-name` (`<run-name>.<stage>.done`), so different phenotype runs sharing
  one `outdir` neither clobber nor wrongly skip each other.
- The sampled-chunk seed is `--seed`; the tool version, the full resolved
  configuration and selected tool versions are recorded in
  `.akgwas/run_manifest[_<run-name>].json`. Keep that file with your results.

## Nothing was written / dry run

`--dry-run` prints every external command it *would* run and changes nothing.
Use it to sanity-check paths before a long build.
