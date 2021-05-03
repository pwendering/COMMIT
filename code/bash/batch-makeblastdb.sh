#!/bin/bash

#########################################################################################################
# Runs makeblastdb (BLAST+, NCBI) on all .faa files in the specified directory                         	#
# Input:                                                                                                #
#               input directory containing the protein fasta files                                      #
#		output directory									#	
# Output:                                                                                               #
#               outputs protein blast databases			                                        #
#########################################################################################################

if [[ "$#" < 2 ]]; then
        echo "USAGE: batch-makeblastdb.sh inputDir outputDir"
        exit 1
fi

# collect all genbank files in an array
input=($(ls $1|grep "\.faa"))

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
	makeblastdb -in $inFile -dbtype prot -title $ID -out $outFile > /dev/null

        echo -e "\t==> done"    
done

