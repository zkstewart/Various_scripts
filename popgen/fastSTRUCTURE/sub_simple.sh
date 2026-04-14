#!/bin/bash -l
#PBS -N fstructure
#PBS -l walltime=24:00:00
#PBS -l mem=50G
#PBS -l ncpus=1
#PBS -J 1-15

cd $PBS_O_WORKDIR

conda activate faststructure

####

# Specify the PLINK2 BED file prefix
## There should be three files like: ${BEDPREFIX}.bed, ${BEDPREFIX}.bim, ${BEDPREFIX}.fam
BEDPREFIX=/mnt/f/lab_members/devindee/reformatted_variants/flgenomics_commercial

# Specify the prefix for output files
PREFIX=flgenomics_commercial

####

SEED=100 # no need to change, a constant seed makes the analysis replicable

# STEP 1: Run structure.py with various K values
structure.py -K ${PBS_ARRAY_INDEX} \
    --input=${BEDPREFIX} --format=bed \
    --output=${PREFIX}_simple \
    --full --seed=${SEED}
