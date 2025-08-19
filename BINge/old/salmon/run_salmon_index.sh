#!/bin/bash -l
#PBS -N salmon_index
#PBS -l walltime=03:30:00
#PBS -l mem=65G
#PBS -l ncpus=8

cd $PBS_O_WORKDIR

## MANUAL SETUP BELOW

# >> SETUP: Specify salmon executable location
SALMONDIR=/home/stewarz2/various_programs/salmon/salmon-1.9.0_linux_x86_64/bin

# >> SETUP: Specify FASTA location
FASTADIR=/home/stewarz2/banana_group/qcav_dge/transcriptome/results/binge_reference
FASTAFILE=qcav_dge_reference.fasta

# >> SETUP: Specify computational resources
CPUS=8

## MANUAL SETUP END

## RUN PROGRAM
# STEP 1: Index the FASTA file where it resides
${SALMONDIR}/salmon index --threads ${CPUS} \
    --transcripts ${FASTADIR}/${FASTAFILE} \
    --index ${FASTAFILE}_salmon_index
