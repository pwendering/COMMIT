#!/bin/bash

#########################################################################################################
# Runs blastp (BLAST+, NCBI) for all .fasta input files against the protein database that corresponds 	#
# to the ID.									                        #
# Input:                                                                                                #
#               input directory containing the protein fasta files					#
#		directory containing the blast protein database						#
#		output directory									#
# Output:                                                                                               #
#               outputs blast results for each BLAST run                                                #
#########################################################################################################

if [[ "$#" < 2 ]]; then
        echo "USAGE: batch-blastp.sh inputDir dbDir outputDir"
        exit 1
fi

# collect all genbank files in an array
input=($(ls $1|grep "\.fasta"))

# create the output directory of it does not exist yet
if [[ ! -d $3 ]]; then
        echo "Creating output directory"
        mkdir -p $3
fi

# Parameters
evalue=10
matrix=BLOSUM90
ncpus=4

# create log file
logFile=$3/batch-blastp.log
echo -e "evalue=$evalue; scoring matrix=$matrix\n" | tee $logFile
echo -e "ID\tOne-many-mappings" >> $logFile

count=0;
for i in ${input[@]};
do
        let count++
        # create the output file name
        ID=$(echo $i | grep -o "^[^\.]*")
        echo "Processing $ID ($count)..."
	inFile=$(echo "$1/$i")
	dbName=$(echo $2/$ID.protDB)
        outFile=$(echo "$3/$ID.aln")
	#echo "$inFile"
	#echo "$outFile"
	#echo "$dbName"
	
	# blastp
	echo -ne "\t> blastp"
	t_1=$(date +%s)
	blastp -db $dbName -query $inFile -out $outFile -evalue $evalue -matrix $matrix -outfmt "6 qseqid sseqid evalue pident" -max_target_seqs 1 -num_threads $ncpus
	echo -e " ($(date -ud "@$(($(date +%s) - $t_1))" +%T))"

	# post process and obtain mapping file
	echo -e "\t> creating tab-separated mapping file"
	#awk -F"\t" -v OFS="\t" 'match($1,/CDS\.[0-9]*/) {print substr($1,RSTART,RLENGTH), $2}' $outFile | sed "s/^/$ID./" | sort | uniq > $3/$ID.mapping
	awk -F"\t" -v OFS="\t" 'match($1,/LOCUS\_[0-9]*/) {print substr($1,RSTART,RLENGTH), $2}' $outFile > $3/$ID.mapping
	# count the amount of duplicate mappings
	n_target=$(cat $3/$ID.mapping | wc -l)
	n_query=$(cat $3/$ID.mapping | cut -f2 | sort | uniq | wc -l)
	echo -e "\t\t $(($n_target-$n_query)) one-to-many mappings"
	echo -e "$ID\t$(($n_target-$n_query))" >> $logFile
        echo -e "\t==> done"
done

