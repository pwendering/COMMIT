#!/bin/bash
#PBS -N "ComGapFill"
#PBS -l nodes=1:ppn=8
#PBS -e ~/ComGapFill/err.txt
#PBS -o ~/ComGapFill/out.txt
matlab -nodesktop -nosplash -r "run_iterative_gap_filling"

