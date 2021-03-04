#!/bin/bash
#PBS -N "ComGapFill_Leaf"
#PBS -l nodes=1:ppn=8
#PBS -e /stud/wendering/Masterthesis/DATA/Gap-filling/iterative/Leaf/all/error200621.txt
#PBS -o /stud/wendering/Masterthesis/DATA/Gap-filling/iterative/Leaf/all/output200621.txt
matlab -nodesktop -nosplash -r "run_iterative_gap_filling"

