#!/bin/bash
#########################################################################################################
# Runs bedtools getfasta on all .gff files in the specified directory                                   #
# Input:                                                                                                #
#               input directory containing the .gff files                                               #
#               input directory containing the .fna files for gff mapping                               #
#               output directory for multi-fasta files                                                  #
# Output:                                                                                               #
#               outputs fasta files from the .gff  files                                                #
#########################################################################################################

if [[ "$#" < 3 ]]
then
        echo "USAGE: batch-convert-gff-to-fasta.sh inputDir fastaDir outputDir"
        exit 1
fi


# collect all gff file names in an array
input=($(ls $1*.gff))

# create the output directory if it does not exist yet
if [[ ! -d $3 ]]
then
        echo "Creating output directory"
        mkdir -p $3
fi

count=0
for i in ${input[@]}
do
        let count++

        # create the output file name
        ID=$(echo $i | grep -o "^[^\.]*")
        echo "Processing $ID ($count) ..."
        echo -e "\t> preparing gff file..."
        inFile=$(echo "$1/$i")
        fastaFile=$(echo "$2/$ID.fna.gz")
        outFile=$(echo "$3/$ID.fasta")
        tmp=$(echo "tmp_$ID")

	# shorten ID in column 9
        awk -F"\t" -v OFS="\t" 'match($9, /CDS.[0-9]+/) {if ($3 ~ /CDS/) {$3=substr($9, RSTART, RLENGTH); print}}' $inFile > $tmp

        # convert gff to fasta
        echo -e "\t> converting..."
        zcat $fastaFile > tmp_fna_$ID
        bedtools getfasta -fi tmp_fna_$ID -bed $tmp -name > $outFile
        echo -e "\t==> done"

        rm $tmp tmp_fna_$ID*
done

