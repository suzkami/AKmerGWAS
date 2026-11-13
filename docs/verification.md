# Verification checklist (run on the Linux box)

This tree was written on Windows with no ability to execute it. Everything
below is meant to be run **once** on the machine that will build/test the conda
package. Work top to bottom; each step is independent enough that a failure
points at a specific layer.

## 0. Environment

```bash
conda create -n akgwas -c conda-forge -c bioconda \
    python=3.12 kmtricks gemma plink bowtie bwa pytest
conda activate akgwas
cd v6
```

## 1. Static checks (no external tools needed)

```bash
python -m compileall akgwas tests          # every module must byte-compile
python -c "import akgwas, akgwas.stages"   # import graph is clean
```

Then the unit tests:

```bash
pytest tests/test_units.py -q
```

This exercises the pure-logic pieces: size parsing, bed pack/unpack round-trip,
bed concat/subset, phenotype reorder/export + ID-mismatch, MAF filter +
chunking, SAM filtering, fof parsing, config defaults + validation.

Watch specifically for:

- `test_bed_concat_and_subset` — confirms our hand-rolled `.bed` byte layout
  matches what PLINK reads. If this fails, the `python` backend is wrong.
- `test_pheno_reorder_and_export` — confirms the ID→fam reordering.

## 2. CLI surface

```bash
python -m akgwas --version          # -> akgwas 1.0.0
python -m akgwas --errors           # prints the E1xx..E9xx table
python -m akgwas build --help       # grouped flags render
python -m akgwas all --help
python -m akgwas all                # must fail: --programhome is required
```

## 3. Synthetic end-to-end

```bash
python -m tests.make_synthetic /tmp/syn
# builds a bowtie index and runs the whole pipeline, then re-runs with --force
pytest tests/test_smoke.py -q
```

The smoke test asserts: the planted `traitA` k-mer appears in
`sigk/traitA_sigK.txt`, the `siggeno`/`samfile`/`sigk_mapped` outputs exist, and
a forced re-run is **byte-identical** (reproducibility). If bowtie is missing
the test is skipped — check it actually ran (`-rs`).

Three fixes land here. The dataset plants a wide panel of **common** k-mers over
20 samples (not 6 samples / 2 variants, and not a panel of rare k-mers either),
so GEMMA's LMM converges instead of dying in GSL `brent` on a rank-deficient or
poorly-conditioned GRM. Only the two **causal** k-mers are kept sparse (at most
four carriers) with **distinct flanking bases**, so no flanking k-mer of a causal
k-mer survives `recurrence-min` to correlate with the trait (previously such
flankers both polluted the GRM and pushed the real causal k-mer out of `sigK`);
flankers of the common panel k-mers are trait-independent and harmless. And the
forced re-run re-enters the `build` stage with `kmatrixs_out/` already present —
that path is the kmtricks run-dir bug fixed in this revision (the stage now
clears and lets kmtricks recreate the directory).

If the run fails, the reason is printed to stderr (with the GEMMA log tail for
`E601`), not swallowed.

## 4. Packaging

```bash
# from the repo root that will be tagged
python -m pip install build && python -m build      # sdist + wheel
python -m pip install dist/akgwas-1.0.0-py3-none-any.whl
akgwas --version
```

For the conda recipe:

```bash
# after tagging v1.0.0 and filling recipe/meta.yaml's sha256:
conda-build v6/recipe/meta.yaml -c conda-forge -c bioconda
```

Notes to confirm while doing this:

- `readme = "README.md"` (capital R) resolves — the 0.1.0 build failed on a
  lowercase `Readme.md`.
- No data files ship in the wheel — the package is pure Python and reads
  nothing from `akgwas/` at runtime (there is no config template in `packages=`).
- The recipe pulls **bowtie** (new default aligner) and drops
  snakemake/perl/rush/parallel.

## 5. Interface checks only a human can do

- **External K/PCA path.** Pass `--use-external-qk --kinship K --pca PCA` with
  a real kinship + PCA you trust, and confirm reordering behaves. Try both
  without `--pca-has-id` (row order assumed) and with `--pca-has-id` (IDs
  reordered/verified).
- **bwa path.** Add `--aligner bwa`, point `ref_list.txt` at bwa index prefixes,
  re-run just `akgwas alignment --programhome ... --force`.
- **plink multibed backend.** Add `--bed-backend plink` and confirm the chunks
  are identical to the `python` backend (`plink --bfile matrix_0001 --recode`
  on both, diff).
- **Real tolerance.** Run on a small real dataset and sanity-check that
  significant k-mers map to plausible reference positions.

## Known review points (please scrutinise)

1. **GEMMA `.assoc.txt` column order** — `kgwas.py` assumes the standard 12
   columns and pulls indices `[1,6,7,8,11]` (rs, af, beta, se, p_wald). Confirm
   against the installed GEMMA version.
2. **gemma `-gk` output name** — `qk.py` now accepts only the real matrix names:
   `QK_kinship.cXX.txt` for `-gk 1` (centred, default) and `QK_kinship.sXX.txt`
   for `-gk 2` (standardized). It never picks up `QK_kinship.log.txt`. Confirm
   your GEMMA build uses exactly those suffixes.
3. **bowtie full-length-exact filter** — `-a -v 0 -f -S` plus the CIGAR `== kM`
   and `NM:i:0` checks in `alignment.py`. Confirm bowtie emits `NM:i:0` for exact
   hits on your version.
4. **plink `--pca` eigenvec shape** — `qk.py` expects FID IID PC1..PC3 with no
   header. Confirm your plink version.
