import argparse
import subprocess
import shutil
from importlib import resources
import os


def get_workflow_path():
    return str(resources.files("akgwas").joinpath("workflow"))


def init_config():
    target = os.path.abspath("config.yaml")

    if os.path.exists(target):
        print("config.yaml already exists")
        return

    with resources.as_file(
        resources.files("akgwas").joinpath("config", "config.yaml")
    ) as config_path:
        shutil.copy(config_path, target)

    print(f"Config file created: {target}")


def run_snakemake(rule_file, configfile, cores):
    workflow_path = get_workflow_path()
    snakefile = os.path.join(workflow_path, "rules", rule_file)

    cmd = [
        "snakemake",
        "-s", snakefile,
        "--configfile", configfile,
        "--cores", str(cores)
    ]

    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        prog="akgwas",
        description="AKmerGWAS: All input format (shor_read; assembly; genomes) for k-mer GWAS analysis.",
        epilog="Example:\n  akgwas --all --configfile config.yaml --cores 8\n  Author: SUZ (https://github.com/suzkami/Aseembly_GWAS)",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--init-config",
        action="store_true",
        help="Initialize a default config.yaml in current directory"
    )

    parser.add_argument("--all",
        action="store_true",
        help="Run full pipeline: build → kgwas → alignment"
    )

    parser.add_argument("--build",
        action="store_true",
        help="Run k-mer counting, matrix building and matrix split step"
    )

    parser.add_argument("--kgwas",
        action="store_true",
        help="Run k-mer based GWAS with multiple threads"
    )

    parser.add_argument("--alignment",
        action="store_true",
        help="Align significant k-mers to reference genomes as list in ref_index_list"
    ) 

    parser.add_argument("--configfile",
        help="!! configuration file, use --init-config to generate defualt"
    )

    parser.add_argument("--cores",
        type=int,
        default=1,
        help="Number of CPU cores to use (default: 1)"
    )

    args = parser.parse_args()

    #init-config
    if args.init_config:
        init_config()
        return

    #config required
    if not args.configfile:
        parser.error("--configfile is required unless using --init-config")

    if not os.path.exists(args.configfile):
        raise FileNotFoundError(f"Config file not found: {args.configfile}")

    #only one module can be specific
    selected = sum([args.build, args.kgwas, args.alignment, args.all])
    if selected != 1:
        parser.error("Please specify exactly one module: --build / --kgwas / --alignment / --all")

    #flow
    if args.all:
        #build → kgwas → alignment
        run_snakemake("build.smk", args.configfile, args.cores)
        run_snakemake("kgwas.smk", args.configfile, args.cores)
        run_snakemake("alignment.smk", args.configfile, args.cores)

    elif args.build:
        run_snakemake("build.smk", args.configfile, args.cores)

    elif args.kgwas:
        run_snakemake("kgwas.smk", args.configfile, args.cores)

    elif args.alignment:
        run_snakemake("alignment.smk", args.configfile, args.cores)


if __name__ == "__main__":
    main()