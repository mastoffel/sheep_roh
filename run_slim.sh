#!/bin/bash

#$ -V
#$ -cwd
#$ -N slim_1000200
#$ -o o_files/
#$ -e e_files/


SCRATCH=/scratch/$USER/$JOB_ID/slim1000200
mkdir -p $SCRATCH

Rscript scripts/slim_sims_pipeline.R $SCRATCH 1000 200

rsync -av $SCRATCH slim_sim/

rm -rf /scratch/$USER/$JOB_ID

