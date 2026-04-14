#!/bin/bash -l
#PBS -N fstructure
#PBS -l walltime=24:00:00
#PBS -l mem=20G
#PBS -l ncpus=1
#PBS -J 1-15

cd $PBS_O_WORKDIR

conda activate faststructure

####

# Specify the location of the PLINK2 BED files
BEDDIR=/home/stewarz2/citrus/devindee/reformatted_variants

# Specify the prefix for input and output files
PREFIX=flgenomics_commercial

## Note: there should be three files within ${BEDDIR} which look like: ${PREFIX}.bed, ${PREFIX}.bim, ${PREFIX}.fam
####

# STEP 1: Run structure.py with various K values
SEED=100 # no need to change, a constant seed makes the analysis replicable
structure.py -K ${PBS_ARRAY_INDEX} \
    --input=${BEDDIR}/${PREFIX} --format=bed \
    --output=${PREFIX}_simple \
    --full --seed=${SEED}
