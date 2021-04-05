#!/bin/bash

#$ -V
#$ -cwd
#$ -N slim_1000200bot
#$ -o o_files/
#$ -e e_files/


SCRATCH=/scratch/$USER/$JOB_ID/slim1000200bot7030del_varyh
mkdir -p $SCRATCH

Rscript scripts/slim_sims_pipeline2.R $SCRATCH 1000 200

rsync -av $SCRATCH slim_sim/

rm -rf /scratch/$USER/$JOB_ID

