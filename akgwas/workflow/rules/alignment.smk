############################################
# alignment.smk
############################################
import os
import glob
from importlib import resources

SCRIPT_DIR = str(resources.files("akgwas").joinpath("workflow", "scripts"))

############################################
WORKDIR = config["programhome"]
PHENO_DIR = os.path.join(WORKDIR, config["pheno_dir"])
OUT_FINAL = os.path.join(WORKDIR, config["output_final"])

############################################
TRAITS = [
    os.path.splitext(os.path.basename(p))[0]
    for p in glob.glob(f"{PHENO_DIR}/*.txt")
]

REFLIST = config["alignment"]["refindexlist"]
REF_DICT = {}

with open(REFLIST) as f:
    for line in f:
        path = line.strip()
        if not path:
            continue
        basename = os.path.basename(path)
        ref_prefix = basename
        for ext in [".fasta.gz.sa", ".fa.gz.sa", ".fasta.sa", ".fa.sa", ".gz.sa", ".sa"]:
            if ref_prefix.endswith(ext):
                ref_prefix = ref_prefix.replace(ext, "")
        REF_DICT[ref_prefix] = path
REFINDEX = list(REF_DICT.keys())

############################################
rule all:
    input:
        expand(
            os.path.join(OUT_FINAL, "samfile/{trait}_{ref_prefix}_matches.sam"),
            trait=TRAITS,
            ref_prefix=REFINDEX
        ),
        expand(
            os.path.join(OUT_FINAL, "sigk_mapped/{trait}_{ref_prefix}_chrpos_sigk.txt"),
            trait=TRAITS,
            ref_prefix=REFINDEX
        ),

############################################
rule gwas_alignment:
    input:
        ref = lambda wc: REF_DICT[wc.ref_prefix],
        pheno = lambda wc: f"{PHENO_DIR}/{wc.trait}.txt"
    output:
        sam = os.path.join(OUT_FINAL, "samfile/{trait}_{ref_prefix}_matches.sam"),
        samchr = os.path.join(OUT_FINAL, "sigk_mapped/{trait}_{ref_prefix}_chrpos_sigk.txt"),
    params:
        work_path = OUT_FINAL,
        kmer_len = config["kmer_len"],
        bwamodel = config["alignment"]["bwamodel"],
    shell:
        """
        bash {SCRIPT_DIR}/kmer_alignment_Large_V3.sh \
            -ref {input.ref} \
            -p {input.pheno} \
            -work_path {params.work_path} \
            -l {params.kmer_len} \
            -map {params.bwamodel}
        """