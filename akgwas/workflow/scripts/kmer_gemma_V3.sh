#!/bin/bash
#limited the threads used to 1 for 1 run
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

set -euo pipefail
# gemma use -lmm 1 ward test
# Required parameter (no default value)
input_bed_prefix=""
input_phenotype=""
input_kinship=""
outtemp_dir="./"
PCA=""
threshold=0.00001

#Parsing command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -bed) input_bed_prefix="$2"; shift ;;
        -p) input_phenotype="$2"; shift ;;
        -k) input_kinship="$2"; shift ;;
        -c) PCA="$2"; shift ;;
        -o) outtemp_dir="$2"; shift ;;
        -threshold) threshold="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

mkdir -p $outtemp_dir

{
if [[ -z "$input_bed_prefix" ]]; then
    echo "Error: must specify input bed path and prefix via -bed"
    exit 1
fi
if [[ -z "$input_phenotype" ]]; then
    echo "Error: must specify input phenotype via -p"
    exit 1
fi
if [[ ! -f "$input_phenotype" ]]; then
    echo "Error: Input phenotype $input_phenotype does not exist"
    exit 1
fi
if [[ -z "$input_kinship" ]]; then
    echo "Error: must specify input kinship via -k"
    exit 1
fi
if [[ ! -f "$input_kinship" ]]; then
    echo "Error: Input kinship $input_kinship does not exist"
    exit 1
fi
if [[ -z "$outtemp_dir" ]]; then
    echo "Error: must specify output file via -o"
    exit 1
fi

cd $outtemp_dir/
prefix=$(basename $input_phenotype | sed 's/.txt//g')
echo -e "gwas for $prefix strat - $outtemp_dir/; start time $(date)"
} > >(tee -a $outtemp_dir/kmer_gemma_aligment.log)

bed_prefix="${input_bed_prefix%.bed}"

#Gemma
gemma -bfile $bed_prefix \
      -k $input_kinship \
      -p $input_phenotype \
      -c $PCA \
      -lmm 1 \
      -o $prefix
{
echo -e "gwas for $prefix done  - Output in $outtemp_dir/output/; end time $(date)"


#combine p and other information and output
echo -e "combine gwas results and other information for $prefix\n  - $outtemp_dir/sigk/"
echo -e "pvalue threhold: $threshold"
mkdir -p $outtemp_dir/sigk

awk -v th="$threshold" 'NR > 1 && $12 <= th' $outtemp_dir/output/${prefix}.assoc.txt |\
    cut -f 2,7-9,12 |\
    sort -k1,1 > $outtemp_dir/sigk/${prefix}_sigK.txt

#prepare sigk fasta, and extrcat bed by plink
awk 'NR > 1 {print $1}' $outtemp_dir/sigk/${prefix}_sigK.txt \
    > $outtemp_dir/sigk/${prefix}_sigrs.txt

mkdir -p $outtemp_dir/siggeno
plink --bfile $bed_prefix \
      --extract $outtemp_dir/sigk/${prefix}_sigrs.txt \
      --make-bed \
      --out $outtemp_dir/siggeno/${prefix}_bed

mkdir -p $outtemp_dir/fasta
awk '{print ">"$1 "\n" $1}' $outtemp_dir/sigk/${prefix}_sigrs.txt \
    > $outtemp_dir/fasta/${prefix}_fasta.txt
rm $outtemp_dir/sigk/${prefix}_sigrs.txt

echo -e "sigk, fasta prepared for $prefix  - Output in $outtemp_dir/sigk/ and $outtemp_dir/fasta/; end time $(date)"
#next step: combime geno; sigk_chrpos_sigk.txt in Main.sh
#then aligment for all phenotypes

} > >(tee -a $outtemp_dir/kmer_gemma_aligment.log)

