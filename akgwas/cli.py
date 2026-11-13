"""Command line interface.

AKmerGWAS is driven entirely from the command line -- there is no config file.

    akgwas build     --programhome /data/run --cores 8
    akgwas kgwas     --programhome /data/run --cores 32
    akgwas alignment --programhome /data/run --cores 8
    akgwas all       --programhome /data/run --cores 32

Every value has a default; pass the matching flag to override it. Only
``--programhome`` is required (it is the base for resolving relative paths and
where outputs are written). ``--run-name`` isolates the phenotype-dependent
outputs so several GWAS runs can share one ``programhome``.
"""

from __future__ import annotations

import argparse
import sys

from . import __version__
from .config import config_from_args
from .errors import AkgwasError, HINTS
from .pipeline import run_module


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="akgwas",
        description="AKmerGWAS: k-mer GWAS for any input format (short reads, assemblies, genomes).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  akgwas build     --programhome /data/run --cores 8\n"
            "  akgwas kgwas     --programhome /data/run --cores 32 --run-name ft10\n"
            "  akgwas alignment --programhome /data/run --aligner bowtie\n"
            "  akgwas all       --programhome /data/run --cores 32\n"
        ),
    )
    parser.add_argument("--version", action="version", version=f"akgwas {__version__}")
    parser.add_argument("--errors", action="store_true",
                        help="print the error-code table and exit")

    sub = parser.add_subparsers(dest="command")

    def add_common(p: argparse.ArgumentParser) -> None:
        p.add_argument("--programhome", required=True,
                       help="REQUIRED. Absolute base directory for paths and outputs.")
        p.add_argument("--cores", type=int, default=1, help="parallel jobs (default: 1)")
        p.add_argument("--dry-run", action="store_true", help="print what would run, change nothing")
        p.add_argument("--force", action="store_true", help="ignore completion markers and re-run stages")
        p.add_argument("-q", "--quiet", action="store_true", help="less console output")

        c = p.add_argument_group("input / output")
        c.add_argument("--outdir", help="output root (default: programhome); relative to programhome")
        c.add_argument("--run-name", help="isolate phenotype-dependent outputs under outdir/<name>/")
        c.add_argument("--fof", help="sample list file, relative to programhome (default: fof.txt)")
        c.add_argument("--pheno", help="phenotype table, relative to programhome (default: pheno.txt)")
        c.add_argument("--kmer-len", type=int, help="k-mer length (default: 31)")
        c.add_argument("--ulimit-min", type=int, help="minimum open-file limit (default: 4096)")
        c.add_argument("--bed-dir", help="shared bed chunk dir name (default: bed)")
        c.add_argument("--qk-dir", help="shared K/PCA dir name (default: QK)")
        c.add_argument("--temp-dir", help="per-run GEMMA temp dir name (default: temp_gemma)")
        c.add_argument("--final-dir", help="per-run result dir name (default: kmergwas_outputs)")
        c.add_argument("--kmtricks-run-dir", help="kmtricks run dir name (default: kmatrixs_out)")

        b = p.add_argument_group("build")
        b.add_argument("--hard-min", type=int, help="min within-sample count to call present (default: 2)")
        b.add_argument("--recurrence-min", type=int, help="cross-sample recurrence min (default: 2)")
        b.add_argument("--partitions", type=int, help="kmtricks partition count; 0 = auto (default: 0)")
        b.add_argument("--build-threads", type=int, help="kmtricks threads (default: 8)")
        b.add_argument("--aggregate-to-stdout", action=argparse.BooleanOptionalAction, default=None,
                        help="stream the aggregate matrix to stdout (default: true)")

        m = p.add_argument_group("multibed")
        m.add_argument("--block-size", help="target chunk size, e.g. 1G (default: 1G)")
        m.add_argument("--min-block-size", help="lower bound for block-size (default: 16M)")
        m.add_argument("--max-block-size", help="upper bound for block-size (default: 10G)")
        m.add_argument("--maf", type=float, help="minor-allele-frequency filter (default: 0.01)")
        m.add_argument("--multibed-threads", type=int, help="multibed threads (default: 8)")
        m.add_argument("--bed-backend", choices=["python", "plink"], help="chunk writer backend (default: python)")

        q = p.add_argument_group("qk (kinship / PCA)")
        q.add_argument("--use-external-qk", action="store_true", default=None,
                       help="use external kinship + PCA instead of the k-mer version")
        q.add_argument("--kinship", help="external kinship matrix (relative to programhome)")
        q.add_argument("--pca", help="external PCA file (relative to programhome)")
        q.add_argument("--sample-fraction", type=float, help="fraction of chunks for internal K/PCA (default: 0.2)")
        q.add_argument("--seed", type=int, help="chunk-sampling seed (default: 123456789)")
        q.add_argument("--gk", type=int, choices=[1, 2], help="gemma -gk: 1 centered, 2 standardized (default: 1)")
        q.add_argument("--pca-has-id", action="store_true", default=None,
                       help="external PCA carries sample IDs in its first columns")
        q.add_argument("--qk-threads", type=int, help="plink --pca threads (default: 1)")

        k = p.add_argument_group("kgwas")
        k.add_argument("--threshold", type=float, help="p_wald extraction threshold (default: 0.0001)")

        a = p.add_argument_group("alignment")
        a.add_argument("--ref-list", help="reference index list, relative to programhome (default: ref_list.txt)")
        a.add_argument("--aligner", choices=["bowtie", "bwa"], help="aligner (default: bowtie)")
        a.add_argument("--max-hits", type=int, help="drop a k-mer with more hits than this (default: 5000)")
        a.add_argument("--alignment-threads", type=int, help="aligner threads (default: 1)")

        co = p.add_argument_group("collect")
        co.add_argument("--gemma-combine", action=argparse.BooleanOptionalAction, default=None,
                        help="also write the full all-k-mers association dump "
                             "(gemma_combine/, can be very large; default: off)")

    for name, help_text in (
        ("build", "k-mer counting (kmtricks) + bed chunking (disk-light streaming)"),
        ("kgwas", "kinship/PCA + per-chunk GEMMA association + result collection"),
        ("alignment", "map significant k-mers to reference genomes (outer loop = reference)"),
        ("all", "run build -> kgwas -> alignment"),
    ):
        sp = sub.add_parser(name, help=help_text, description=help_text)
        add_common(sp)

    return parser


def print_errors() -> None:
    for code in sorted(HINTS):
        print(f"{code}  {HINTS[code]}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.errors:
        print_errors()
        return 0

    if not args.command:
        parser.print_help()
        return 0

    try:
        cfg = config_from_args(args)
        return run_module(
            cfg,
            args.command,
            cores=args.cores,
            dry_run=args.dry_run,
            force=args.force,
            verbose=not args.quiet,
        )
    except AkgwasError as exc:
        print(exc.render(), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
