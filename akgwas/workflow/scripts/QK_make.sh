#!/bin/bash
#build kinship and pca from bed/ file

bedinpath="/home/suz/kmer_home/snaketest/bed"
bedoutpath="/home/suz/kmer_home/snaketest/QK"
NumberbedSelect=0.2

#Parsing command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -bedpath) bedinpath="$2"; shift ;;
        -o) bedoutpath="$2"; shift ;;
        -bfb) NumberbedSelect="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

mkdir -p $bedoutpath
Numberbeds=$(ls "$bedinpath"/*.bed 2>/dev/null | wc -l)

#random extract NumberbedSelect*Numberbeds for kinship and pca calculate
to_select=$(awk -v n="$Numberbeds" -v s="$NumberbedSelect" 'BEGIN{print int(n*s + (n*s > int(n*s) ? 1 : 0))}')
if [ "$to_select" -lt 1 ]; then
    to_select=1
fi
mapfile -t random_numbers < <(shuf -i 1-"$Numberbeds" -n "$to_select" | sort -n)

mergelist="$bedoutpath/mergelist.txt"
echo -n "" > "$mergelist" 

for num in "${random_numbers[@]}"; do
    ls $bedinpath/matrix_${num}.{bed,bim,fam} | xargs >> "$mergelist"
done

#PCA (Q)
plink --merge-list $mergelist --make-bed --out $bedoutpath/QK_bed
plink --bfile $bedoutpath/QK_bed --pca 3 --out $bedoutpath/QK_pca
awk '{print $3, $4, $5}' $bedoutpath/QK_pca.eigenvec > $bedoutpath/QK_pca1-3.txt

#Kinship (K)
#give false pheno to makesure gemma accept correct format
awk '{$6=1; print}' "$bedoutpath/QK_bed.fam" > "$bedoutpath/QK_bed.fam.tmp" && \
mv "$bedoutpath/QK_bed.fam.tmp" "$bedoutpath/QK_bed.fam"

cd $bedoutpath/
gemma -bfile $bedoutpath/QK_bed \
      -gk -miss 1.0 -maf 0.0 -r2 1.0 \
      -o QK_kinship

mv $bedoutpath/output/* $bedoutpath/
rm -rf $bedoutpath/output/
rm $bedoutpath/QK_bed.*

#OUTPUT in $bedoutpath/