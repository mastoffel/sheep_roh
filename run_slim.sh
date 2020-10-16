#!/bin/bash

#$ -V
#$ -cwd
#$ -N slim_sim_test
#$ -o o_files/
#$ -e e_files/

SCRATCH=/scratch/$USER/$JOB_ID/slim200
mkdir -p $SCRATCH

Rscript scripts/slim_combined.R $SCRATCH 200 200

rsync -av $SCRATCH ./

rm -rf /scratch/$USER/$JOB_ID

