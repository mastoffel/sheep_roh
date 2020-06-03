#!/bin/sh	
# Grid engine options
#$ -N transpose_ped
#$ -wd /exports/csce/eddie/biology/groups/pemberton/martin/sheep_roh/
#$ -l h_rt=03:00:00
#$ -M martin.adam.stoffel@gmail.com
#$ -m beas
#$ -l h_vmem=60G
#$ -o o_files/
#$ -e e_files/


/exports/csce/eddie/biology/groups/pemberton/martin/utils/bin/garlic/bin/linux/garlic --tped data/sheep_geno_trans.tped --tfam data/sheep_geno_trans.tfam --centromere data/dummy_cent --error 0.005 --winsize 60 --out output/roh 
