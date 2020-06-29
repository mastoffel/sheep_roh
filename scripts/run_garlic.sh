#!/bin/sh	
# Grid engine options
#$ -N run_garlic
#$ -wd /exports/csce/eddie/biology/groups/pemberton/martin/sheep_roh/
#$ -l h_rt=00:30:00
#$ -M martin.adam.stoffel@gmail.com
#$ -m beas
#$ -l h_vmem=30G
#$ -o o_files/
#$ -e e_files/


/exports/csce/eddie/biology/groups/pemberton/martin/utils/bin/garlic/bin/linux/garlic --tped data/sheep_geno_filt.tped --tfam data/sheep_geno_filt.tfam --centromere data/dummy_cent --error 0.005 --winsize 100 --overlap-frac 0.5 --out output/roh 
