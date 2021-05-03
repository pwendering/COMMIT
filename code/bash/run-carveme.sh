#!/bin/bash
# Reconstruct genome-scale metabolic model(s) from genomic sequences using CarveMe

if [[ $# != 1 ]]; then
  echo -e "Usage: run-carveme type [possible types: 'prot', 'nucl']"
fi

# First check if carveMe is installed:
installCheck=$(which carve 2>&1)
if [[ ! -z $(echo $installCheck | grep "no carve") ]];
then
echo -e "\nCarveMe is either not installed, not on the PATH or you do not have the permission to use it.\n"
fi

# run CarveMe on all peptide or DNA FASTA files in the current directory

if [[ $1 == "prot" ]];
then
  opt=""
  type="faa"
elif [[ $1 == "nucl" ]];
then
  opt="--dna"
  type = "fna"
fi


while read line
do
  name=$(echo $line | sed "s/.$type//")
  echo -e "Processing $name...\n\n"
  # carving
  echo -e "Carving process for $line...\n"
  carve -r $opt $line -o $name.xml.gz
done < <(ls | grep $(echo "*.$type") )



