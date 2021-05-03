#!/bin/bash
#################################################################################
# Runs hmmbuild for all files in the given directory containing 		#
# multiple sequence alignments (MSAs)						#
# Input:									#
# 		input directory	containing ONLY the MSA files			#
#		output directory						#
# Output:									#
#		outputs the Hidden Markov models in the specified output folder	#
#################################################################################

if [[ "$#" < 2 ]]; then
	echo "USAGE: run-batch-hmmbuild.sh inputDir outputDir"
	exit 1
fi

# save input file names in an array
input=($(ls $1))

# create the output directory of it does not exist yet
if [[ ! -d $2 ]]; then
	echo "Creating output directory"
	mkdir $2
fi

for i in ${input[@]};
do	
	# create the output file name
	ID=$(echo $i | grep -o "^[^\.]*")
	echo $ID
	outFile=$(echo "$ID.hmm")
	nonempty=$(cat $1/$i | grep -m 1 "[a-z]")
	if [[ -f "$outFile" ]]; then
		echo -e "\t==> HMM for $ID already exists"
	elif [[ -z $nonempty ]]; then
		echo -e "\t==> no MSA"
	else
		# run hmmbuild
		hmmbuild -o out.txt $2/$outFile $1/$i
		
	fi
done
