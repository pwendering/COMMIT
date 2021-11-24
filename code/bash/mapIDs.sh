#!/bin/bash
####################################################################### 
### Search for an exact match of a search pattern and extract 
### the desired column mathing to this key
### Input:  $1: file containing one column containing the key (source ids)
###		    $2: column of the target namespace to be extracted
###		    $3: file name of the translation file to be used
### Output: List of translated ids (or "") on stdout
####################################################################### 

# Column containing the target namespace
col_val=$2

while read line; do
	# Escape special characters before using grep 
	query=$(echo $line | sed -e 's/+/\\+/g' -e 's/-/\\-/g')
	# Before and after the pattern are either white space or '\|'
	match=$(grep -E "(^|\||[[:space:]])+$query([[:space:]]|\|)+" $3 | cut -f$col_val)
  if [[ -z $match ]]
  then
    echo $query
  else
    echo $match
  fi
    
done < $1 	
	


