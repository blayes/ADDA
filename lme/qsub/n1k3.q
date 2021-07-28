#!/bin/bash
#$ -N n1k3
#$ -q BLAYES
#$ -pe smp 24
#$ -wd /Users/ssrivastva/dist_lme/code/
#$ -t 1-40
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load stack/2020.1
module load r-matrix/1.2-17_gcc-9.2.0 r-matrixstats/0.55.0_gcc-9.2.0 

R CMD BATCH --no-save --no-restore "--args 1 3 $SGE_TASK_ID" submit_part.R samp/n1k3_$SGE_TASK_ID.rout


