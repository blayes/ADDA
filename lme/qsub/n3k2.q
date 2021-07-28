#!/bin/bash
#$ -N n3k2
#$ -q INFORMATICS
#$ -pe smp 4
#$ -wd /Users/ssrivastva/dist_lme/code/
#$ -t 1-40
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load stack/2020.1
module load r-matrix/1.2-17_gcc-9.2.0 r-matrixstats/0.55.0_gcc-9.2.0 

R CMD BATCH --no-save --no-restore "--args 3 2 $SGE_TASK_ID" submit_part.R samp/n3k2_$SGE_TASK_ID.rout


