#!/bin/bash -l
#PBS -N salmon_index
#PBS -l walltime=04:00:00
#PBS -l mem=50G
#PBS -l ncpus=4

cd $PBS_O_WORKDIR

## MANUAL SETUP BELOW

# >> SETUP: Specify salmon executable location
SALMONDIR=/home/stewarz2/various_programs/salmon/salmon-1.9.0_linux_x86_64/bin

# >> SETUP: Specify FASTA location
FASTADIR=/home/stewarz2/anemones/cassie/transcriptome/final_transcriptome
FASTAFILE=cassie_okay-okalt.fasta

# >> SETUP: Specify computational resources
CPUS=1

## MANUAL SETUP END

## RUN PROGRAM
# STEP 1: Index the FASTA file where it resides
${SALMONDIR}/salmon index --threads ${CPUS} \
	--transcripts ${FASTADIR}/${FASTAFILE} \
	--index ${FASTAFILE}_salmon_index
