#!/bin/bash -l
#PBS -N cactus
#PBS -l walltime=32:00:00
#PBS -l mem=120G
#PBS -l ncpus=16

cd $PBS_O_WORKDIR

conda activate cactus

####

# Specify parameters
CPUS=16
PREFIX=14Q021_18Q053_pangenome

####

JOBSTOREPATH=${PWD}/cactus_intermediate
SEQFILE=${PWD}/cactus_samples.txt
OUTDIR=${PWD}/cactus_output
REFERENCE=reference

# STEP 1: Run cactus-pangenome
cactus-pangenome ${JOBSTOREPATH} ${SEQFILE} \
    --outDir ${OUTDIR} --outName ${PREFIX} \
    --reference ${REFERENCE} --maxCores ${CPUS} \
    --binariesMode singularity \
    --vcf --vcfwave --giraffe
