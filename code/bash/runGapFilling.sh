#!/bin/bash
#PBS -N "ComGapFill"
#PBS -l nodes=1:ppn=8
#PBS -M philipp.wendering@gmail.com
#PBS -m ae
#PBS -e /stud/wendering/ComGapFill/err.txt
#PBS -o /stud/wendering/ComGapFill/out.txt
matlab -nodesktop -nosplash -r "gap_fill_individual_models"

