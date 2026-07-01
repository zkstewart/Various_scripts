#!/bin/bash -l
#PBS -N cactus
#PBS -l walltime=48:00:00
#PBS -l mem=120G
#PBS -l ncpus=16

cd $PBS_O_WORKDIR

conda activate cactus

#################################

# Specify parameters
CPUS=316
PREFIX=glauca # output directory name
REFERENCE=14Q021.1 # sample as specified in 'cactus_samples.txt' that should be treated as the reference/backbone of the pangenome

#################################

JOBSTOREPATH=${PWD}/cactus_intermediate
SEQFILE=${PWD}/cactus_samples.txt
OUTDIR=${PWD}/cactus_output

# STEP 1: Run cactus-pangenome
cactus-pangenome ${JOBSTOREPATH} ${SEQFILE} \
    --outDir ${OUTDIR} --outName ${PREFIX} \
    --reference ${REFERENCE} --maxCores ${CPUS} \
    --binariesMode singularity \
    --vcf --vcfwave --giraffe
