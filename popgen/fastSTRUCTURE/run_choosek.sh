#!/bin/bash -l
#PBS -N choosek
#PBS -l walltime=00:30:00
#PBS -l mem=10G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

conda activate faststructure

####

# Specify the location of the fastSTRUCTURE outputs (with various K values)
RESULTSDIR=/mnt/f/lab_members/devindee/structure

# Specify the prefix for the fastSTRUCTURE output files
PREFIX=flgenomics_commercial

####

# STEP 1: Run script to identify optimal K for explaining population structure
chooseK.py --input=${RESULTSDIR}/${PREFIX}
