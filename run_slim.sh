#!/bin/bash

#$ -V
#$ -cwd
#$ -N slim_const_200
#$ -o o_files/
#$ -e e_files/

SCRATCH=/scratch/$USER/$JOB_ID/slim200const
mkdir -p $SCRATCH

Rscript scripts/slim_combined.R $SCRATCH 200 200

rsync -av $SCRATCH slim_sim/

rm -rf /scratch/$USER/$JOB_ID

