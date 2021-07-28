#!/bin/bash
#$ -N n1k2
#$ -q all.q
#$ -pe smp 8
#$ -wd /Users/ssrivastva/dist_mov_len/code/
#$ -t 1-40
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load stack/2021.1
module load r-matrix/1.3-2_gcc-9.3.0
module load r-matrixstats/0.58.0_gcc-9.3.0

R CMD BATCH --no-save --no-restore "--args 1 2 $SGE_TASK_ID" submit_part_logistic_eps.R samp/n1k2_$SGE_TASK_ID.rout






