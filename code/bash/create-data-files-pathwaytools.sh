#!/bin/bash
########################################################################################################
# Write all required attributes for a list of samples that are required for the generation of PGDBs 
# by PathoLogic in batch mode and copies the assembly and annotatinos files to the respective folders
# Input: 
#		Folder containing assemblies
#		Folder containing the annotations to the assemblies
# 		File containing the taxonomic units of the the samples (left column sample names that
#			correspond to the names of the assembly and annotation files, right column 
#			taxonomic units)
#
# Output:
#		folders containing the required inputs for PathoLogic
########################################################################################################
if [[ $# < 2 ]]; then
	echo "USAGE: assemblyDir annotationDir taxonomyDir"	
fi

taxFile=$3

annot_suffix=".gbk"
assbl_suffix=".fna"

arrayIDs=($(ls -R $1 | grep $assbl_suffix))

echo -e "\n\n"
echo "Assembly directory: $1"
echo "Annotation directory: $2"
echo "Taxonomy file: $taxFile"
echo -e "\n\n"

for i in "${arrayIDs[@]}"
do
	# get all required values
	ID=$(echo $i | grep -o "^[^\.]*")
	NAME=$(grep -o "^$ID\s" $taxFile| cut -f2)
	STORAGE="File"
	AUTHOR="Philipp Wendering:Potsdam University"
	DBNAME=$(echo "${ID}DB")
	
	echo -e "Sample:\t$ID\t($NAME)"

	# create the organism-params file
        dname=$(echo $ID)
        fname=$(echo "$dname/organism-params.dat")
        mkdir $dname
        touch $fname

	# write the values to the new file
	echo "#$ID data file for generation of PGDB files with PathoLogic" >> $fname
	echo -e "ID\t$ID" >> $fname
	echo -e "STORAGE\t$STORAGE" >> $fname
	echo -e "NAME\t$NAME" >> $fname
	echo -e "AUTHOR\t$AUTHOR" >> $fname
	echo -e "DBNAME\t$DBNAME" >> $fname
	
	# create the genetic-elements file
	fname=$(echo "$dname/genetic-elements.dat")
	echo -e "ID\t$ID" >> $fname
	echo -e "ANNOT-FILE\t$ID$annot_suffix" >> $fname
	echo -e "SEQ-FILE\t$ID$assbl_suffix" >> $fname

	# copy the assembly file uncompress if neccessary
	habitat=$(echo $ID | grep -o "^[^0-9]*")
	cp $1/$habitat/$i $dname
	if [[ ! -z $(echo $i | grep "\.gz") ]]; then
		gunzip $dname/$i
	fi
	
	# copy the annotation file to the organisms' folder and add the organism name
	cp $2/$habitat/$ID$annot_suffix $dname/$ID$annot_suffix
	sed -i -e "s/ORGANISM/ORGANISM  $NAME/g" -e "s/SOURCE/SOURCE      $NAME./g" "$dname/$ID$annot_suffix"

	# create script.lisp
	touch $dname/script.lisp
	orgID=$(echo $ID | tr [a-z] [A-Z])
	echo "(in-package :ecocyc)" >> $dname/script.lisp
	echo "(select-organism :org-id '$orgID)" >> $dname/script.lisp
	echo "(create-flat-files-for-current-kb)" >> $dname/script.lisp
	
done

# create file to store all directory names for Pathway tools batch mode
ls -d $PWD/* > param-dirs-list.txt
