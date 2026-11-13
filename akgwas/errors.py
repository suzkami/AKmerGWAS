"""Centralised error codes for AKmerGWAS.

Every user-facing failure funnels through :class:`AkgwasError`. The error code
(``E###``) is stable and is meant to be looked up in ``docs/troubleshooting.md``
so that a bare message is never the only clue a user gets.

The codes are grouped by area:

    E1xx  environment / preflight (tools, ulimit, paths)
    E2xx  input validation (fof, phenotype, reference list, config)
    E3xx  sample alignment / sample-order invariants
    E4xx  build + multibed (kmtricks, chunking, bed writing)
    E5xx  kinship / PCA (plink, gemma -gk)
    E6xx  GWAS (gemma -lmm)
    E7xx  collection / merge
    E8xx  alignment (bowtie / bwa / SAM parsing)
    E9xx  internal / unexpected
"""

from __future__ import annotations

# --- E1xx environment ------------------------------------------------------
E_TOOL_MISSING = "E101"
E_ULIMIT_TOO_LOW = "E102"
E_WORKDIR_UNWRITABLE = "E103"

# --- E2xx input validation -------------------------------------------------
E_INPUT_MISSING = "E201"
E_CONFIG_INVALID = "E202"
E_FOF_MALFORMED = "E204"
E_ROW_COUNT_MISMATCH = "E205"
E_REF_LIST_EMPTY = "E206"
E_PHENO_MALFORMED = "E207"

# --- E3xx sample alignment -------------------------------------------------
E_SAMPLE_ID_MISMATCH = "E301"
E_SAMPLE_ORDER_MISMATCH = "E302"
E_SAMPLE_ORDER_UNVERIFIABLE = "E303"

# --- E4xx build / multibed -------------------------------------------------
E_KMTricks_FAILED = "E401"
E_EMPTY_MATRIX = "E402"
E_BED_BACKEND_UNKNOWN = "E403"
E_MATRIX_MALFORMED = "E404"

# --- E5xx kinship / PCA ----------------------------------------------------
E_PLINK_FAILED = "E501"
E_GEMMA_GK_FAILED = "E502"

# --- E6xx GWAS -------------------------------------------------------------
E_GEMMA_FAILED = "E601"
E_ASSOC_MALFORMED = "E602"

# --- E7xx collect ----------------------------------------------------------
E_COLLECT_EMPTY = "E701"

# --- E8xx alignment --------------------------------------------------------
E_ALIGNER_UNKNOWN = "E801"
E_ALIGN_FAILED = "E802"

# --- E9xx internal ---------------------------------------------------------
E_INTERNAL = "E999"

#: Human readable hints, keyed by code. Printed alongside the error message.
HINTS: dict[str, str] = {
    E_TOOL_MISSING: "Install the missing tool (conda env) and make sure it is on PATH.",
    E_ULIMIT_TOO_LOW: "Raise the open-file limit, e.g. `ulimit -n 4096` before running.",
    E_WORKDIR_UNWRITABLE: "Check directory permissions and free disk space.",
    E_INPUT_MISSING: "Check the flag/path; relative inputs resolve against --programhome.",
    E_CONFIG_INVALID: "Check the value; run `akgwas <subcommand> --help` for valid ranges and defaults.",
    E_FOF_MALFORMED: "fof.txt format is `sample : file1 ; file2 ; ...` (one sample per line).",
    E_ROW_COUNT_MISMATCH: "The phenotype file must contain one row per sample in fof.txt.",
    E_SAMPLE_ID_MISMATCH: "Phenotype/fof sample IDs must form the same set as the bed .fam IDs.",
    E_SAMPLE_ORDER_MISMATCH: "Kinship/PCA must be built with the same sample order as the GWAS bed.",
    E_EMPTY_MATRIX: "kmtricks produced no k-mers; check --hard-min / --recurrence-min and the inputs.",
    E_MATRIX_MALFORMED: "A matrix row did not match fof sample count or contained non-0/1 values.",
    E_PLINK_FAILED: "Inspect the captured plink log; often a malformed bed/bim/fam.",
    E_GEMMA_FAILED: "Inspect the captured gemma log; usually a kinship/phenotype format issue.",
    E_ASSOC_MALFORMED: "The GEMMA .assoc.txt had an unexpected column count.",
    E_ALIGNER_UNKNOWN: "`--aligner` must be one of: bowtie, bwa.",
    E_INTERNAL: "This is a bug; please report it with the run manifest.",
}


class AkgwasError(Exception):
    """A user-facing pipeline error tagged with a stable code."""

    def __init__(self, code: str, message: str, *, hint: str | None = None) -> None:
        self.code = code
        self.message = message
        self.hint = hint if hint is not None else HINTS.get(code)
        super().__init__(message)

    def render(self) -> str:
        parts = [f"[{self.code}] {self.message}"]
        if self.hint:
            parts.append(f"    hint: {self.hint}")
        return "\n".join(parts)
