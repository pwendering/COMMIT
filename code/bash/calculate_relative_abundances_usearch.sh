#!/bin/bash
# Use usearch to calculate realtive abundances of OTUs given by a fasta file of 16S rRNA sequences
# Inputs:
# 	$1	sample to analyze
#	$2 	16S rRNA seqences
#	$3 	identity level [0,1]
#	$4 	output directory

set -e

if [[ $# < 4 ]]; then
	echo "USAGE: calculate_relative_abundances_usearch sample-to-analyze 16S-rRNA-sequences identity-level output-directory"
	exit 1
fi

# rename input arguments
sample=$1
OTUseqs=$2
identity=$3
outDir=$(echo $4 | sed 's/\/$//')

# dereplicate sample
if [[ ! -f $outDir/seqs_unique.fasta ]]; then
	echo "-----------------------------------------------------------"
	echo " Dereplicating sample..."
	echo "-----------------------------------------------------------"
	usearch -fastx_uniques $sample \
		-fastaout $outDir/seqs_unique.fasta \
		-sizeout
fi

# sort by abundance and remove singletons
#echo -e "\ncalculating abundances with a similarity value of $(($identity * 100))%...\n"
#usearch -closed_ref $outDir/seqs_unique.fasta \
#	-db $OTUseqs \
#	-id $identity \
#	-strand both \
#	-tabbedout $outDir/closed_ref_$(($identity * 100)).txt

## extract numbers for the single OTUs
#echo -e "\nextacting abundances...\n"
#cat $outDir/closed_ref_$(($identity * 100)).txt | cut -f 4 | sort | uniq -c | sort -n | grep -v "\*" > $outDir/abundances_sorted.tsv

# create OTU table
echo "-----------------------------------------------------------"
echo " Creating OTU table..."
echo "-----------------------------------------------------------"
usearch -otutab $outDir/seqs_unique.fasta \
       	-otus $OTUseqs \
       	-id $identity \
       	-otutabout $outDir/otutab.txt \
       	-mapout $outDir/map.txt

# normalize OTU table and remove rare OTUs
#echo "-----------------------------------------------------------"
#echo " Normalizing OTU table with usearch..."
#echo "-----------------------------------------------------------"
#usearch -otutab_rare $outDir/otutab.txt \
#       	-sample_size 500 \
#      	-output $outDir/otutab_rare.txt

# converting to biom format
echo "-----------------------------------------------------------"
echo " Converting to biom format... "
echo "-----------------------------------------------------------"
biom convert -i $outDir/otutab.txt \
	-o $outDir/otutab.biom \
	--table-type="OTU table" \
	--to-json

# CSS normalization
echo "-----------------------------------------------------------"
echo " Performing CSS normalization..."
echo "-----------------------------------------------------------"
$HOME/ComGapFill/CSS_normalization.R \
	$outDir/otutab.biom \
	$outDir/otutab_norm.biom 2> /dev/null

# convert normalized table to text format
echo "-----------------------------------------------------------"
echo " Converting normalized table to text format"
echo "-----------------------------------------------------------"
biom convert -i $outDir/otutab_norm.biom \
	-o $outDir/otutab_norm.txt \
	--table-type="OTU table" \
	--to-tsv
echo "-----------------------------------------------------------"
echo "done"
