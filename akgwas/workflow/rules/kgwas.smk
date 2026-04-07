############################################
# kgwas.smk
############################################
import os
import glob
from importlib import resources

SCRIPT_DIR = str(resources.files("akgwas").joinpath("workflow", "scripts"))

############################################
WORKDIR = config["programhome"]
BED_CHUNK = os.path.join(WORKDIR, "bed")
PHENO_DIR = os.path.join(WORKDIR, config["pheno_dir"])
QK_OUT = os.path.join(WORKDIR, config["QK_path"])
GEMMA_TEMP = os.path.join(WORKDIR, config["kwags_temp"])
OUT_FINAL = os.path.join(WORKDIR, config["output_final"])

############################################
TRAITS = [
    os.path.splitext(os.path.basename(p))[0]
    for p in glob.glob(f"{PHENO_DIR}/*.txt")
]

CHUNKS = [
    os.path.splitext(os.path.basename(b))[0]
    for b in glob.glob(f"{BED_CHUNK}/matrix_*.bed")
]

############################################
rule all:
    input:
        #outpath: (OUT_FINAL, 1-5)
        # 1. sigk
        expand(
            os.path.join(OUT_FINAL, "sigk/{trait}_sigK.txt"),
            trait=TRAITS
        ),
        # 2. fasta
        expand(
            os.path.join(OUT_FINAL, "fasta/{trait}_fasta.txt"),
            trait=TRAITS
        ),
        # 3. bed
        expand(
            os.path.join(OUT_FINAL, "siggeno/{trait}_bed.bed"),
            trait=TRAITS
        ),
        # 4. gemma.assoc
        expand(
            os.path.join(OUT_FINAL, "gemma_combine/{trait}.assoc_combined.txt"),
            trait=TRAITS
        ),
        # 5. gemma.log
        expand(
            os.path.join(OUT_FINAL, "logs/{trait}.log_combined.txt"),
            trait=TRAITS
        )


############################################
USE_EXTERNAL_QK = config["kgwas"].get("use_external_qk", False)
USER_KINSHIP = config["kgwas"]["kinship"]
USER_PCA = config["kgwas"]["pca"]

if USE_EXTERNAL_QK:
    KINSHIP = USER_KINSHIP
    PCA = USER_PCA
    QK_CMD = f"""
        echo "Using external Q/K matrix files"
        echo "Kinship: {USER_KINSHIP}"
        echo "PCA: {USER_PCA}"
        
        # 验证文件存在
        if [ ! -f "{USER_KINSHIP}" ]; then
            echo "ERROR: Kinship file not found: {USER_KINSHIP}"
            exit 1
        fi
        if [ ! -f "{USER_PCA}" ]; then
            echo "ERROR: PCA file not found: {USER_PCA}"
            exit 1
        fi
        
        ln -sf {USER_KINSHIP} {KINSHIP}
        ln -sf {USER_PCA} {PCA}
    """
else:
    KINSHIP = os.path.join(QK_OUT, "QK_kinship.cXX.txt")
    PCA = os.path.join(QK_OUT, "QK_pca1-3.txt")
    QK_CMD = """
        echo "Running qkmatrix..."
        bash {SCRIPT_DIR}/QK_make.sh \
            -bedpath {input} \
            -o {params.outdir} \
            -bfb {params.select}
    """

############################################
rule qkmatrix:
    input:
        BED_CHUNK
    output:
        kinship = KINSHIP,
        pca = PCA
    params:
        outdir = QK_OUT,
        select = config["kgwas"]["NumberbedSelect"]
    shell:
        QK_CMD

############################################
rule gwas_chunk:
    input:
        bed = lambda wildcards: os.path.join(BED_CHUNK, f"{wildcards.chunk}.bed"),
        pheno = lambda wildcards: os.path.join(PHENO_DIR, f"{wildcards.trait}.txt"),
        kinship = KINSHIP,
        pca = PCA
    output:
        sigk = os.path.join(GEMMA_TEMP, "{chunk}/sigk/{trait}_sigK.txt"),
        fasta = os.path.join(GEMMA_TEMP, "{chunk}/fasta/{trait}_fasta.txt"),
        assoc = os.path.join(GEMMA_TEMP, "{chunk}/output/{trait}.assoc.txt"),
        log = os.path.join(GEMMA_TEMP, "{chunk}/output/{trait}.log.txt")
    threads: config["kgwas"].get("threads", 1)
    params:
        threshold = config["kgwas"]["threshold"],
    shell:
        """
        bash {SCRIPT_DIR}/kmer_gemma_V3.sh \
            -bed {input.bed} \
            -p {input.pheno} \
            -k {input.kinship} \
            -c {input.pca} \
            -threshold {params.threshold} \
            -o {GEMMA_TEMP}/{wildcards.chunk}
        """

############################################
rule kgwas_sum:
    input:
        lambda wc: expand(
            rules.gwas_chunk.output.sigk,
            chunk=glob_wildcards(f"{BED_CHUNK}/{{chunk}}.bed").chunk,
            trait=wc.trait
        )
    output:
        sigk = os.path.join(OUT_FINAL, "sigk/{trait}_sigK.txt"),
        fasta = os.path.join(OUT_FINAL, "fasta/{trait}_fasta.txt"),
        assoc = os.path.join(OUT_FINAL, "gemma_combine/{trait}.assoc_combined.txt"),
        log = os.path.join(OUT_FINAL, "logs/{trait}.log_combined.txt"),
        bed = os.path.join(OUT_FINAL, "siggeno/{trait}_bed.bed")
    params:
        outtemp = GEMMA_TEMP,
        outdir = OUT_FINAL,
        phenodir = PHENO_DIR
    shell:
        """
        bash {SCRIPT_DIR}/gwasresul_collecet.sh \
            -temp {params.outtemp} \
            -o {params.outdir} \
            -p {params.phenodir}
        """