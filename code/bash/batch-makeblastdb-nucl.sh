#!/bin/bash

#########################################################################################################
# Runs makeblastdb (BLAST+, NCBI) on all .fna files in the specified directory                         	#
# Input:                                                                                                #
#               input directory containing the nucleotide fasta files                                   #
#		output directory									#	
# Output:                                                                                               #
#               outputs nucleotide blast databases			                                #
#########################################################################################################

if [[ "$#" < 2 ]]; then
        echo "USAGE: batch-makeblastdb-nucl.sh inputDir outputDir"
        exit 1
fi

# collect all multi-fasta nucleotide sequence files in an array
input=($(ls $1|grep "\.fna"))

# create the output directory of it does not exist yet
if [[ ! -d $2 ]]; then
        echo "Creating output directory"
        mkdir $2
fi

count=0;
for i in ${input[@]};
do
        let count++
        # create the output file name
        ID=$(echo $i | grep -o "^[^\.]*")
        echo "Processing $ID ($count)..."
        inFile=$(echo "$1/$i")
        outFile=$(echo "$2/$ID.protDB")
	#echo "$inFile"
	#echo "$outFile"
        # makeblastdb
	makeblastdb -in $inFile -dbtype nucl -title $ID -out $outFile > /dev/null

        echo -e "\t==> done"    
done

