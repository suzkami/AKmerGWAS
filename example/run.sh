#a pipline run the full stream of akgwas
PATHWORK=$(pwd)

mkdir $PATHWORK/sequence
mkdir $PATHWORK/bwaindex

wget -i arabi_sequence_url.txt -c -w 5 -P $PATHWORK/sequence
find $PATHWORK/sequence -type f -name "*.fasta.gz" | head -n 3 | rush -j 3 "bwa index {} -p $PATHWORK/bwaindex/{%}"
find $PATHWORK/bwaindex -type f -name "*.sa" > ref_list.txt

#because config.yaml already prepared by author, here skip
#akgwas --init-config
#vim config.yaml
sed -i "s|^programhome: \"/home/suz/kmer_home/AKmerGWAS/example\"|programhome: \"$PATHWORK\"|" config.yaml
akgwas --all --configfile config.yaml --cores 8