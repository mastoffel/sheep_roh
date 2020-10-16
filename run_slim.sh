#!/bin/bash

#$ -V
#$ -cwd
#$ -N slim_sim_10000
#$ -o o_files/
#$ -e e_files/

SCRATCH=/scratch/$USER/$JOB_ID/slim10000
mkdir -p $SCRATCH

Rscript scripts/slim_combined.R $SCRATCH 10000 200

rsync -av $SCRATCH slim_sim/

rm -rf /scratch/$USER/$JOB_ID

