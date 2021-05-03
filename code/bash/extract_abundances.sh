#!/bin/bash
# extract the OTU abundances of a clustering file resulting from usearch
# Input
#	$1 input file resulting from usearch pipeline
# Output
#	file containing the abundance of each OTU respectively (absolute)

includeTiess=0
while getopts 't' opt; do
	case $opt in
		t) includeTies=1 ;;
		*) echo "Error in parsing flags"
	esac
done

shift "$(( OPTIND - 1 ))"

if [[ "$includeTies" -eq 1 ]]; then
	cat $1 | cut -f 4,6 | sed -e 's/ties=//g' -e 's/[0-9]\+://g' | grep -oE "\w+" | sort | uniq -c | sort -n | grep -v "\<0\>"
else
	cat $1 | cut -f 4 | sort | uniq -c | sort -n | grep -v "\*"
fi

