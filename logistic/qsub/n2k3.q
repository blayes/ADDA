#!/bin/bash
#$ -N n2k3
#$ -q LT
#$ -pe smp 56
#$ -wd /Users/ssrivastva/dist_logistic/code/
#$ -t 1-40
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load stack/2020.1
module load r-matrix/1.2-17_gcc-9.2.0 r-matrixstats/0.55.0_gcc-9.2.0 

R CMD BATCH --no-save --no-restore "--args 2 3 $SGE_TASK_ID" submit_part.R samp/n2k3_$SGE_TASK_ID.rout


