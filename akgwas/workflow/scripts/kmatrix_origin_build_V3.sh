#!/bin/bash

#requires package kmtricks
#requires input fastq file as fofile list: each sample and fastq one line, if paired, just separate
    #e.g. (HOR_448 : /home/a_R1.fastq.gz ; /home/a_R2.fastq.gz)
set -euo pipefail

#variable
ksize=31        #kmer mer
k_ind_min=2     #lower limit of k-mer count calculated for each sample
k_all_min=2     #lower limit of k-mer frequence calculated for all sample
thread=8        #threads number
nb_partitions=0 #number of partitions (0=auto)

# Required parameter (no default value)
input_file=""
output_dir=""


#Parsing command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -l) ksize="$2"; shift ;;
        -n) k_ind_min="$2"; shift ;;
        -t) thread="$2"; shift ;;
        -fof) input_file="$2"; shift ;;
        -o) output_dir="$2"; shift ;;
        -r) k_all_min="$2"; shift ;;
        -p) nb_partitions="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

{
if [[ -z "$input_file" ]]; then
    echo "Error: must specify input file via -file"
    exit 1
fi
if [[ -z "$output_dir" ]]; then
    echo "Error: must specify output file via -o"
    exit 1
fi
if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file $input_file does not exist"
    exit 1
fi

echo "===== Parameter configuration ====="
echo "kmer_size: $ksize"
echo "hard_min: $k_ind_min"
echo "recurrence_min: $k_all_min"
echo "Threads: $thread"
echo "input_file: $input_file"
echo "output_dir: $output_dir"
echo "Final matrix path: $output_dir/bed"
echo "number of partitions: $nb_partitions"

#step 1: build martix
start_time=$(date +%s)
echo -e "start kmer matrix count: $(date)"
} > >(tee -a $output_dir/kmatrix_build.log)

#seperate count and merge, merge with large threads will case error (I/O problem)
kmtricks pipeline --file $input_file \
                  --run-dir $output_dir/kmatrixs_out \
                  --kmer-size $ksize \
                  --mode kmer:pa:bin \
                  --hard-min $k_ind_min \
                  --recurrence-min $k_all_min \
                  --nb-partitions $nb_partitions \
                  -t $thread \
                  --until count \
                  --cpr

kmtricks merge --run-dir $output_dir/kmatrixs_out \
               --recurrence-min $k_all_min \
               --mode kmer:pa:bin \
               --cpr \
               -t 2
{
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo -e "kmer matrix count done, end time: $(date); time used: ${elapsed_time}" 

######
#step 2: merge and dumo martix
mkdir -p "$output_dir/Matrix"

start_time=$(date +%s)
echo -e "start kmer pa matrix build: $(date)"
kmtricks aggregate --run-dir $output_dir/kmatrixs_out --pa-matrix kmer --cpr-in -t $thread \
        --format text --output $output_dir/Matrix/Origin_matrix.txt

#when Origin_matrix.txt constructed, delet $output_dir/kmatrixs_out
rm -rf $output_dir/kmatrixs_out

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo -e "kmer pa matrix build done, end time: $(date)\ntime used: ${elapsed_time}"
} > >(tee -a $output_dir/kmatrix_build.log)