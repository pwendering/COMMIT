#!/bin/bash
#########################################################################################################
# Run PathoLogic on a list of given input folders							#
# Input: 												#
#		File containing the absolute paths to the folders					#
#													#
# Output:												#
#		.dat files in the respective input files that are required for AuReMe reconstruction	#
#		output file that summarizes the progress						#
#########################################################################################################

if [[ $# < 1 ]]; then
	echo "USAGE: run-pathway-tools dirFile"
	exit 1
elif [[ ! -f "$1" ]]; then
	echo "File not found $1"
	exit 1
fi

pgdbDir=pathway/ptools-local/pgdbs/user

# create output file
touch ptools-out.txt

# create the PGDB files
echo "#################### Creating PGDB files ####################" >> ptools-out.txt
while read line;
do	
	ID=$(echo $line | grep -o "[^/]*$")
	echo -e "Sample:\t$ID" >> ptools-out.txt
	ID=$(echo $ID | tr '[:upper:]' '[:lower:]')
	ID=$(echo "${ID}cyc")

	if [[ ! -d $pgdbDir/$ID ]]; then
		pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -patho $line
		echo -e "\t==> Done." >> ptools-out.txt
	else
		echo -e "\t==> PGDB already exists." >> ptools-out.txt
	fi

done < $1

# create .dat files
echo "#################### Creating .dat files ####################" >> ptools-out.txt
while read line;
do
	ID=$(echo $line | grep -o "[^/]*$")
	ID=$(echo $ID | tr '[:upper:]' '[:lower:]')
	ID=$(echo "${ID}cyc")
	echo -e "Sample:\t$ID" >> ptools-out.txt
	if [[ -f $pgdbDir/$ID/1.0/data/proteins.dat ]]; then
		echo -e "\t==> .dat files already exist" >> ptools-out.txt
	elif [[ -d $pgdbDir/$ID ]]; then
		pathway-tools -no-web-cel-overview -no-cel-overview -no-patch-download -disable-metadata-saving -nologfile -load $line/script.lisp &
		sleep 15s
		PT=$!
		c=0
		while true;
		do
			n=$(ls $pgdbDir/$ID/1.0/data | grep dat | wc -l)
			if [ $n == "23" ] || [c == "10"]; then
				kill $PT
				break
			else
				sleep 30s
				echo "waiting 30s"
				let "c++"
			fi
		done

		echo -e "\t==> Done." >> ptools-out.txt
	else
		echo -e "\t==> no PGDB found" >> ptools-out.txt
	fi
	
done < $1

# copy .dat files
echo "#################### Copying .dat files ####################" >> ptools-out.txt
while read line;
do
	ID=$(echo $line | grep -o "[^/]*$")
	ID=$(echo $ID | tr '[:upper:]' '[:lower:]')
	ID=$(echo "${ID}cyc")
	if [[ -d $pgdbDir/$ID ]]; then
		echo -e "Sample:\t$ID" >> ptools-out.txt
		if [[ ! -f $line/classes.dat ]]; then
			fileDir=$pgdbDir/$ID/1.0/data
			targetDir=$line
			if [[ ! -d $targetDir ]]; then 
				mkdir $targetDir
			fi
			cp $fileDir/proteins.dat $targetDir
			cp $fileDir/reactions.dat $targetDir
			cp $fileDir/genes.dat $targetDir
			cp $fileDir/enzrxns.dat $targetDir
			cp $fileDir/pathways.dat $targetDir
			cp $fileDir/compounds.dat $targetDir
			cp $fileDir/classes.dat $targetDir
			cp $fileDir/metabolic-reactions.xml $targetDir
			echo -e "\t==> Done." >> ptools-out.txt
		else
			echo -e "\t==> .dat files already exist." >> ptools-out.txt
		fi
	else
		echo -e "\t==> no PGDB found" >> ptools-out.txt
	fi

done < $1


