#!/bin/bash

set -euo pipefail

thread=8
blocksize=10G
fof_file=""
input_file=""
output_dir=""

#Parsing command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -fof) fof_file="$2"; shift ;;
        -file) input_file="$2"; shift ;;
        -o) output_dir="$2"; shift ;;
        -m) maf="$2"; shift ;;
        -b) blocksize="$2"; shift ;;
        -t) thread="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

{
start_time=$(date +%s)
echo "Start maf filter: start time: $(date)"

#Run here
mkdir -p $output_dir/temp
mkdir -p $output_dir/bed

Num_samples=$(wc -l < $fof_file)
limit_count=$(awk -v n=$Num_samples -v maf=$maf 'BEGIN{print int(n - maf*n)}')
cut -d ' ' -f1 $fof_file > $output_dir/bed/sample_id_inorder.txt

parallel --pipe -j $thread --block $blocksize --max-lines 0 '
awk -v limit='"$limit_count"' -v OFS=" " "
{
    zero_count=0
    one_count=0
    for (i=2;i<=NF;i++) {
        if (\$i==0) zero_count++
        else if (\$i==1) one_count++
    }
    if (zero_count<=limit && one_count<=limit) {
        gsub(\"1\",\"2 2\")
        gsub(\"0\",\"1 1\")
        print \"1\",\$1,\"0\",\"0\",\$0
    }
}
" | cut -d "' '" -f1-4,6- > '$output_dir'/temp/matrix_{#}.tped
' < "$input_file"

cut -d ' ' -f 1 $fof_file > $output_dir/bed/sample_id_inorder.txt

#tped tfam
for F in $(find "$output_dir/temp/" -name "matrix_*.tped" -exec basename {} .tped \;)
do
    awk '{print "0 "$1" 0 0 0 0" }' $output_dir/bed/sample_id_inorder.txt \
        > $output_dir/temp/${F}.tfam
done

#transfer each piece to bed
find "$output_dir/temp/" -name "matrix_*.tped" -exec basename {} .tped \; |\
    rush -j $thread "plink --tfile $output_dir/temp/{} --make-bed --silent --out $output_dir/bed/{}"

find "$output_dir/bed/" -name "matrix_*.bed" -exec basename {} .bed \; |\
    rush -j $thread "echo 'After filter, total rows of {} files: ' && wc -l \"$output_dir/bed/{}.bim\" | awk '{print \$1}'"

#when matrix_*.bed constructed, delet $output_dir/temp/
rm -rf $output_dir/temp/

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo -e "bed build build done, end time: $(date)\ntime used: ${elapsed_time}"

} > >(tee -a $output_dir/kmatrix_to_bed.log)