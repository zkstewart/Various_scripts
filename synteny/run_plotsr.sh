#!/bin/bash -l
#PBS -N plotsrPipe
#PBS -l walltime=00:15:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

conda activate synteny

####

# Specify the Various_scripts location
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the location of the input files
GENDIR=/home/stewarz2/plant_group/juel/glauca_progeny/synteny/citrus/chr2_citrus_sequences

GEN1=${GENDIR}/reticulata_chr2.fasta
GEN2=${GENDIR}/sinensis_1_chr2.fasta
GEN3=${GENDIR}/limon_chr2.fasta

# Specify how many CPUs to use
CPUS=16

# Specify prefix for outputs
PREFIX=chr2_ft3_plotsr

####

# STEP 1: Run the pipeline
python ${VARSCRIPTDIR}/synteny/plotsr_pipeline.py --threads ${CPUS} \
    -i ${GEN1} ${GEN2} ${GEN3} \
    -o ${PREFIX} \
    --width 7 --height 10

#--chr <chr1 chr2>
#--markers </location/of/markers.bed>
#--reg <chr2:start-end>
