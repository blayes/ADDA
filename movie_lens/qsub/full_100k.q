#!/bin/bash
#$ -N full_100k
#$ -q LT
#$ -pe smp 4
#$ -wd /Users/ssrivastva/dist_mov_len/code/
#$ -t 1-10
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load stack/2021.1
module load r-matrix/1.3-2_gcc-9.3.0

R CMD BATCH --no-save --no-restore "--args 3 $SGE_TASK_ID" ml_submit_full.R full/ml1_$SGE_TASK_ID.rout
