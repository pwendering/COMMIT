#!/bin/bash
#########################################################################################################
# Run AuReMe (padmet utils: https://gitlab.inria.fr/maite/padmet-utils) on a list of PGDB sets		#
# to create SBML V3 files										#
													#
# Input: 												#
#		File containing the path to the folders							#
#		output top directory									#
#													#
# Output:												#
#		.sbml files containing the generated metabolic models					#
#########################################################################################################

# padmet-utils must be located in $HOME/bin/padmet-utils
UTILPATH=$HOME/bin/padmet-utils/connection
# SBML level
SBML_LVL=3
# MNXref chem_prop directory
MNX_CHEMPROP=$HOME/COMMIT/data/tables/MNXref/chem_prop.tsv
# MNXref chem_xref directory
MNX_CHEMXREF=$HOME/COMMIT/data/tables/MNXref/chem_xref.tsv

if [[ "$#" < 2 ]]; then
	echo "USAGE: run-AuReMe.sh param-dirs-list outTopDir"
	exit 0
elif [[ ! -d $UTILPATH ]]; then
	echo "Please make sure padmet-utils are installed and located in \$HOME/bin/Fedora"
	exit 1
elif [[ ! -f $1 ]]; then
	echo "The given input file does not exist"
elif [[ ! -d $2 ]]; then
	echo "The given output top directory does not exist, do you with to create it? [y|n]"
	read response
	if [[ $response == "y" ]]; then
		mkdir -p $2
	else
		exit 0
	fi
fi

touch $2/error.txt

while read line;
do
	ID=$(echo $line | grep -o "[^/]*$")
	echo -e "Sample:\t$ID"
	habitat=$(echo $ID | grep -oi "^[a-z]*")
	outDir=$2/$habitat
	mkdir $outDir 2> /dev/null
	if [[ -f $outDir/$ID.sbml ]]; then
		echo -e "\t==> sbml file already exists."
	elif [[ -z $(ls $line | grep "classes.dat") ]]; then
		echo -e "\t==> no .dat files found."
	else
		# creating padmet file
		python3.7 $UTILPATH/pgdb_to_padmet.py --version=22.0 --db=$ID --output=$line/$ID.padmet --directory=$line -g -m 1> /dev/null 2>> $2/error.txt
		
		# creating sbml file
		python3.7 $UTILPATH/sbmlGenerator.py --padmet=$line/$ID.padmet --output=$outDir/$ID.sbml --sbml_lvl=$SBML_LVL --mnx_chem_prop=$MNX_CHEMPROP --mnx_chem_xref=$MNX_CHEMXREF 1> /dev/null 2>> $2/error.txt

		echo -e "\t==> Done."
	fi

done < $1

