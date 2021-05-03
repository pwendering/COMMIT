#!/bin/bash
# iterates of all subdirectories that contain a file 'seqs.fasta' (processed reads)
echo -e "\nIterating over the following directories containing \'seqs.fasta\':\n"
find -name "seqs.fasta"
echo ""
sleep 5s

for file in $(find -name "seqs.fasta")
do
	wd=$(echo $file | sed 's/\/seqs.fasta//')
	echo ""
	echo "-----------------------------------------------------------"
	echo "Current working directory: $wd"
	echo ""
	~/COMMIT/code/bash/calculate_relative_abundances_usearch.sh $wd/seqs.fasta $wd/../*isolates*.fasta 1 $wd

done


