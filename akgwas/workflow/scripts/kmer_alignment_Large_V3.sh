#!/bin/bash
#map sig kmer to genomes
#input need:
    #kmer_gemma.sh output fasta file
    #and reference genome
set -euo pipefail
# Required parameter (no default value)
input_refgenome=""
work_path=""
input_phenotype=""
ksize=""
map="bwa"

#Parsing command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -ref) input_refgenome="$2"; shift ;;
        -work_path) work_path="$2"; shift ;;
        -p) input_phenotype="$2"; shift ;;
        -l) ksize="$2"; shift ;;
        -map) map="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

{
if [[ -z "$ksize" ]]; then
    echo "Error: must specify kmer size via -l"
    exit 1
fi

mkdir -p $work_path/samfile
mkdir -p $work_path/sigk_mapped

#mapping
prefix=$(basename $input_phenotype | sed 's/.txt//g')
refindex="${input_refgenome%.sa}"
ref_prefix=$(basename $input_refgenome | sed 's/.fasta.gz.sa//g' | sed 's/.fa.gz.sa//g' | sed 's/.fasta.sa//g' | sed 's/.gz.sa//g' | sed 's/.sa//g')

echo -e "start aligment ${prefix} to ${ref_prefix}; start time $(date)"

${map} mem -k $ksize -T $ksize \
    -a $refindex \
    $work_path/fasta/${prefix}_fasta.txt > $work_path/samfile/${prefix}_${ref_prefix}_matches.sam

sort -k1,1 $work_path/sigk/${prefix}_sigK.txt \
    > $work_path/sigk/${prefix}_sigK_sorted.txt

awk '$2 != 4' $work_path/samfile/${prefix}_${ref_prefix}_matches.sam |\
    cut -f 1,2,3,4,5,6 |\
    sort -k1,1 > $work_path/sigk_mapped/${prefix}_matchedsam_sorted_1-6col.txt

join -1 1 -2 1 -t $'\t' $work_path/sigk_mapped/${prefix}_matchedsam_sorted_1-6col.txt \
    $work_path/sigk/${prefix}_sigK_sorted.txt |\
    sort -k3,3 -k4,4n \
    > $work_path/sigk_mapped/${prefix}_${ref_prefix}_chrpos_sigk.txt

rm $work_path/sigk_mapped/${prefix}_matchedsam_sorted_1-6col.txt
echo -e "${prefix} aligment to ${ref_prefix} done; job done; end time: $(date)"

} > >(tee -a $work_path/kmer_aligment.log)
