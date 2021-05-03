#!/bin/bash

# Create Translation tables from the MetaNetX universal database
dbDir="$HOME/COMMIT/data/tables/MNXref"
metDB=$dbDir/chem_xref.tsv
metpropDB=$dbDir/chem_prop.tsv
rxnDB=$dbDir/reac_xref.tsv
rxnpropDir=$dbDir/reac_prop.tsv
R=$dbDir/MNXref_RXNS.csv
M=$dbDir/MNXref_METS.csv

echo ""
######################### Create two-column Files that contain the matching IDs ######################### 

###### METABOLITES

echo -e "Generating metabolite translation files\n"

# only metabolite keys
cat $dbDir/$metpropDB | grep -v "#" | cut -f1 > $M

# KEGG
echo "MNXref ==> KEGG"
cat $dbDir/$metDB | grep -v "#" | cut -f1,2 | grep "kegg:[CG]" | sed 's/kegg://' | awk -v OFS="\t" '{print $2,$1}' > $dbDir/KEGG_METS_FROM_MNXref.csv

# BiGG
echo "MNXref ==> BiGG"
cat $dbDir/$metDB | grep -v "#" | cut -f1,2 | grep "bigg:" | sed 's/bigg://' | awk -v OFS="\t" '{print $2,$1}' > $dbDir/BiGG_METS_FROM_MNXref.csv

# MetaCyc
echo "MNXref ==> MetaCyc"
cat $dbDir/$metDB | grep -v "#" | cut -f1,2 | grep "metacyc:" | sed 's/metacyc://' | awk -v OFS="\t" '{print $2,$1}' > $dbDir/MetaCyc_METS_FROM_MNXref.csv

# SEED
echo "MNXref ==> SEED"
cat $dbDir/$metDB | grep -v "#" | cut -f1,2 | grep "seed:" | sed 's/seed://' | awk -v OFS="\t" '{print $2,$1}' > $dbDir/ModelSEED_METS_FROM_MNXref.csv

# ChEBI
echo "MNXref ==> ChEBI"
cat $dbDir/$metDB | grep -v "#" | cut -f1,2 | grep "chebi:" | sed 's/chebi://' | awk -v OFS="\t" '{print $2,$1}' > $dbDir/ChEBI_METS_FROM_MNXref.csv

# names
echo "MNXref metabolite names"
cat $dbDir/$metpropDir | grep -v "#" | awk -F "\t" -v OFS="\t" '{print $1,$2}' > $dbDir/MET_NAMES_FROM_MNXref.csv

###### REACTIONS

echo -e "Generating reaction translation files\n\n"
# only reaction keys
cat $dbDir/$rxnpropDir | grep -v "#" |  cut -f1

# KEGG
echo "MNXref ==> KEGG"
cat $dbDir/$rxnDB | grep -v "#" | cut -f1,2 | grep "kegg:" | sed 's/kegg://' | awk -v OFS="\t" '{print $2,$1}' > $dbDir/KEGG_RXNS_FROM_MNXref.csv

# BiGG
echo "MNXref ==> BiGG"
cat $dbDir/$rxnDB | grep -v "#" | cut -f1,2 | grep "bigg:" | sed 's/bigg://' | awk -v OFS="\t" '{print $2,$1}' > $dbDir/BiGG_RXNS_FROM_MNXref.csv

# MetaCyc
echo "MNXref ==> MetaCyc"
cat $dbDir/$rxnDB | grep -v "#" | cut -f1,2 | grep "metacyc:" | sed 's/metacyc://' | awk -v OFS="\t" '{print $2,$1}' > $dbDir/MetaCyc_RXNS_FROM_MNXref.csv

# SEED
echo "MNXref ==> SEED"
cat $dbDir/$rxnDB | grep -v "#" | cut -f1,2 | grep "seed:" | sed 's/seed://' | awk -v OFS="\t" '{print $2,$1}' > $dbDir/ModelSEED_RXNS_FROM_MNXref.csv

# Rhea
echo "MNXref ==> Rhea"
cat $dbDir/$rxnDB | grep -v "#" | cut -f1,2 | grep "rhea:" | sed 's/rhea://' | awk -v OFS="\t" '{print $2,$1}' > $dbDir/Rhea_RXNS_FROM_MNXref.csv

# Enzyme Commission numbers
echo "MNXref ==> E.C."
cat $dbDir/$rxnpropDir | grep -v "#" |  cut -f1,5 > $dbDir/EC_FROM_MNXref.csv


######################### Create translation tables that contain all IDs mapped to MNXref namespace #########################

###### REACTIONS

mapFiles.sh $R $dbDir/KEGG_RXNS_FROM_MNXref.csv > $dbDir/KEGG_col
mapFiles.sh $R $dbDir/BiGG_RXNS_FROM_MNXref.csv > $dbDir/BiGG_col
mapFiles.sh $R $dbDir/MetaCyc_RXNS_FROM_MNXref.csv > $dbDir/MetaCyc_col
mapFiles.sh $R $dbDir/ModelSEED_RXNS_FROM_MNXref.csv > $dbDir/SEED_col
mapFiles.sh $R $dbDir/Rhea_RXNS_FROM_MNXref.csv > $dbDir/Rhea_col
mapFiles.sh $R $dbDir/EC_FROM_MNXref.csv > $dbDir/EC_col

# Combine and write to file
header="MNXref\tKEGG\tBiGG\tMetaCyc\tModelSEED\tRhea\tEC\n"
echo -e $header > $dbDir/MNXref-rxn-translation-table.csv
paste $R $dbDir/KEGG_col BiGG_col $dbDir/MetaCyc_col $dbDir/SEED_col $dbDir/Rhea_col EC_col >> $dbDir/MNXref-rxn-translation-table.csv

####### METABOLITES

mapFiles.sh $M $dbDir/KEGG_METS_FROM_MNXref.csv > $dbDir/KEGG_col
mapFiles.sh $M $dbDir/BiGG_METS_FROM_MNXref.csv > $dbDir/BiGG_col
mapFiles.sh $M $dbDir/MetaCyc_METS_FROM_MNXref.csv > $dbDir/MetaCyc_col
mapFiles.sh $M $dbDir/ModelSEED_METS_FROM_MNXref.csv > $dbDir/SEED_col
mapFiles.sh $M $dbDir/ChEBI_METS_FROM_MNXref.csv > $dbDir/ChEBI_col
mapFiles.sh $M $dbDir/MET_NAMES_FROM_MNXref.csv > $dbDir/NAMES_col

header="MNXref\tKEGG\tBiGG\tMetaCyc\tModelSEED\tChEBI\tNAMES\n"
echo -e $header > $dbDir/MNXref-met-translation-table.csv
paste $M $dbDir/KEGG_col $dbDir/BiGG_col $dbDir/MetaCyc_col $dbDir/SEED_col $dbDir/NAMES_col >> $dbDir/MNXref-met-translation-table.csv

rm $dbDir/KEGG_col $dbDir/BiGG_col $dbDir/MetaCyc_col $dbDir/SEED_col $dbDir/Rhea_col $dbDir/EC_col $dbDir/ChEBI_col $dbDir/NAMES_col

echo -e "\nDONE.\n"
