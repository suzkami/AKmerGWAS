# AKmerGWAS

AKmerGWAS is a pipeline for k-mer based genome-wide association studies (GWAS) using raw sequencing data or assembled genomes.

---

## 📦 Installation

Install via conda:

```bash
conda install -c bioconda akgwas
```

Test if the installation is successful:

```bash
akgwas --help
```

or

```bash
akgwas --version
```

---

## 🚀 Quick Start

```bash
git clone https://github.com/suzkami/AKmerGWAS
cd AKmerGWAS/example
bash run.sh
```

---

## 📂 Input Requirements

To run AKmerGWAS, prepare the following files:

### 1. Sequence Data (FOF file)

Raw sequencing files or assembled genomes. Create a `fof.txt` file (one line per sample).

**Format:**
```
sample_name : file1.fasta.gz ; file2.fasta.gz ; ... ; fileN.fasta.gz
```

**Example:**
```
A88 : /path/to/pseudo88_1.fasta.gz ; /path/to/pseudo88_2.fasta.gz
A108 : /path/to/pseudo108.fasta.gz
```

---

### 2. Phenotype Files

- Each phenotype stored in a separate `.txt` file  
- Single column format with no title 
- Missing values filled with `NA`  
- Must match the sample order in `fof.txt`  

**Example:**
```
61.25
70
68.75
72
63.5
```

---

### 3. Configuration File

Generate template:

```bash
akgwas --init-config
```

Then modify `config.yaml` as needed.

---

### 4. Reference Genome Index (BWA)

Prepare BWA index files for the reference genome.

**Example files:**
```
pseudo10001.fasta.gz.amb
pseudo10001.fasta.gz.ann
pseudo10001.fasta.gz.bwt
pseudo10001.fasta.gz.pac
pseudo10001.fasta.gz.sa
```

Create a `ref_list.txt` file listing `.sa` files:

```
/path/to/pseudo10001.fasta.gz.sa
/path/to/pseudo10002.fasta.gz.sa
```

---

### 5. Optional: Kinship and PCA Files

- Must match sample number and order in `fof.txt`
- If not provided, AKmerGWAS will generate them automatically
- Results may vary due to random sampling

To ensure reproducibility, back up generated files in:
```
QK_path: "QK"
```

---

## ⚠️ Notes

- Recommended population size: **at least 96 samples**
- RAM required: **depends on block_size and threads used, 1GB block_size per thread need approximate 3GB RAM**
- The main bottleneck: **high storage space requirements (for short-read data, at least 3 times the data size is required)**
- Example for 96 Arabidopsis samples: **25 GB of RAM (1GB block_size and 8 threads used) and approximately 150 GB of storage**
- Run for 1135 arabidopsis: **50GB RAM and ~1.5TB storge required**

---

## 📁 Working Directory Structure

```
/home/kmer_home/kmergwas_test
├── config.yaml
├── fof.txt
├── pheno
│   ├── Arabi1135_ft10.txt
│   └── Arabi1135_ft16.txt
└── ref_list.txt
```

---

## Usage

Initialize config:

```bash
akgwas --init-config
```

Edit configuration:

```bash
vim config.yaml
```

Run full pipeline:

```bash
akgwas --all --configfile config.yaml --cores 8
```

Or run steps separately:

```bash
akgwas --build --configfile config.yaml --cores 8
akgwas --kgwas --configfile config.yaml --cores 8
akgwas --alignment --configfile config.yaml --cores 8
```

---

## Input Examples

### Reference Sequences

```
/home/sdc/Reference/Arabi/1001/
├── pseudo10001.fasta.gz
├── pseudo10002.fasta.gz
```

---

### fof.txt

```
head -n 2 /home/kmer_home/kmergwas_test/fof.txt 
A88 : /home/sdc/Reference/Arabi/1001/pseudo88.fasta.gz
A108 : /home/sdc/Reference/Arabi/1001/pseudo108.fasta.gz
```

---

### BWA Index

```
/home/sdc/Reference_index/Arabi/
├── pseudo10001.fasta.gz.amb
├── pseudo10001.fasta.gz.ann
├── pseudo10001.fasta.gz.bwt
├── pseudo10001.fasta.gz.pac
├── pseudo10001.fasta.gz.sa
```

---

### ref_list.txt

```
head -n 2 /home/kmer_home/kmergwas_test/ref_list.txt
/home/sdc/Reference_index/Arabi/pseudo10001.fasta.gz.sa
/home/sdc/Reference_index/Arabi/pseudo10002.fasta.gz.sa
```

---

### Phenotype File

```
head -n 5 /home/kmer_home/kmergwas_test/pheno/Arabi1135_ft10.txt
61.25
70
68.75
72
63.5
```

---

## 📤 Output

### --build

```
Matrix                 # Origin large k-mer Matrix
bed/Matrix_chunks      # multi chunks of Origin large k-mer Matrix, bed format, used for donwstream kgwas
```

---

### --kgwas

```
QK/
kmergwas_outputs/
├── sigk/              # {pheno}_sigk.txt (significant k-mers); Titles: QNAME(k-mer) af  beta  se  p_wald
├── siggeno/           # {pheno}_bed.{bed,bim,fam} (PLINK format)
├── fasta/             # sequences for alignment
├── logs/              # log files
├── gemma_combine/     # GEMMA raw outputs
```

---

### --alignment

```
kmergwas_outputs/
├── samfile/           # sam file, sigk k-mer mapped to each reference
├── sigk_mapped/
│   # {pheno}_{reference}_chrpos_sigk.txt
│   # mapped k-mer results
│   # Titles: QNAME(k-mer)  FLAG  RNAME  POS  MAPQ  CIGAR  af  beta  se  p_wald
```
