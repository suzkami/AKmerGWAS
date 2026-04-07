############################################
# build.smk
############################################

#写死os.path.join(config["programhome"], config["build_matrix"], "Origin_matrix.txt") --> os.path.join(config["programhome"], "Matrix/Origin_matrix.txt")
#写死os.path.join(config["programhome"], config["build_bed_chunks"]) --> os.path.join(config["programhome"], "bed")

import os
import glob
from importlib import resources

SCRIPT_DIR = str(resources.files("akgwas").joinpath("workflow", "scripts"))

############################################
rule all:
    input:
       # matrix = os.path.join(config["programhome"], "Matrix/Origin_matrix.txt"),
       bed_dir = directory(os.path.join(config["programhome"], "bed"))

rule build:
    input:
        fof = config["fof_file"]
    output:
        matrix = os.path.join(config["programhome"], "Matrix/Origin_matrix.txt")
    threads:
        config["build"]["threads"]
    params:
        outdir = config["programhome"],
        kmer_len = config["kmer_len"],
        hard_min = config["build"]["hard_min"],
        partitions = config["build"]["partitions"]
    shell:
        """
        bash {SCRIPT_DIR}/kmatrix_origin_build_V3.sh \
            -l {params.kmer_len} \
            -n {params.hard_min} \
            -t {threads} \
            -p {params.partitions} \
            -fof {input.fof} \
            -o {params.outdir}
        """


#################################################
# multibed.smk
#################################################
rule multibed:
    input:
        matrix = os.path.join(config["programhome"], "Matrix/Origin_matrix.txt")
    output:
        directory(os.path.join(config["programhome"], "bed"))
    threads:
        config["multibed"]["threads"]
    params:
        outdir = config["programhome"],
        fof = config["fof_file"],
        block_size = config["multibed"]["block_size"],
        maf = config["multibed"]["maf"]
    shell:
        """
        bash {SCRIPT_DIR}/matrix_bed_multithreads_V3.sh \
            -fof {params.fof} \
            -file {input.matrix} \
            -o {params.outdir} \
            -b {params.block_size} \
            -m {params.maf} \
            -t {threads}
        """