#!/bin/bash
#build kinship and pca from bed/ file

OUTTEMP="/home/suz/kmer_home/snaketest/kmergwas_outtemps"
PHENOPATH="/home/suz/kmer_home/snaketest/Testpheno"
OUTHOME="/home/suz/kmer_home/snaketest/kmergwas_outputs"

#Parsing command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -temp) OUTTEMP="$2"; shift ;;
        -p) PHENOPATH="$2"; shift ;;
        -o) OUTHOME="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

mkdir -p $OUTHOME/sigk
mkdir -p $OUTHOME/fasta
mkdir -p $OUTHOME/siggeno
mkdir -p $OUTHOME/logs
mkdir -p $OUTHOME/gemma_combine
#mkdir -p $OUTHOME/samfile
#mkdir -p $OUTHOME/sigk_mapped

sleep 180
for P in $(basename -a $PHENOPATH/*.txt);
do  
    Index=$(echo ${P} | sed 's/.txt//g')
    for chunk in $(find $OUTTEMP/ -maxdepth 1 -mindepth 1 -type d -printf "%f\n")
    do
        cat $OUTTEMP/${chunk}/sigk/${Index}_sigK.txt \
        >> $OUTHOME/sigk/${Index}_sigK.txt

        cat $OUTTEMP/${chunk}/fasta/${Index}_fasta.txt \
        >> $OUTHOME/fasta/${Index}_fasta.txt

        echo $OUTTEMP/${chunk} >> $OUTHOME/logs/${Index}.log_combined.txt
        cat $OUTTEMP/${chunk}/output/${Index}.log.txt \
        >> $OUTHOME/logs/${Index}.log_combined.txt

        cat $OUTTEMP/${chunk}/output/${Index}.assoc.txt |\
        cut -f 2,7-9,12 \
        >> $OUTHOME/gemma_combine/${Index}.assoc_combined.txt

        ls $OUTTEMP/${chunk}/siggeno/${Index}_*.{bed,bim,fam} | xargs \
        >> $OUTHOME/siggeno/${Index}_mergelist.txt
    done

    plink --merge-list $OUTHOME/siggeno/${Index}_mergelist.txt --make-bed --out $OUTHOME/siggeno/${Index}_bed
done

#rm -rf $OUTTEMP/