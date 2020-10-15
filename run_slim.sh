#!/bin/bash

#$ -V
#$ -cwd
#$ -N slim_sim
#$ -o o_files/
#$ -e e_files/

SCRATCH=/scratch/$USER/$JOB_ID/slim
mkdir -p $SCRATCH


