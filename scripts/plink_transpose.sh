#!/bin/sh	
# Grid engine options
#$ -N transpose_ped
#$ -wd /exports/csce/eddie/biology/groups/pemberton/martin/sheep_roh/
#$ -l h_rt=00:10:00
#$ -M martin.adam.stoffel@gmail.com
#$ -m beas
#$ -l h_vmem=10G
#$ -pe sharedmem 4
#$ -o o_files/
#$ -e e_files/

# Initialise the environment modules
. /etc/profile.d/modules.sh

# load plink
module load igmm/apps/plink/1.90b4

plink --bfile data/sheep_geno_imputed_oar31_17052020 --sheep --out data/sheep_geno_trans --recode --transpose --exclude data/oar_imp_low_call_snps095.txt
