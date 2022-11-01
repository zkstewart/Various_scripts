#!/bin/bash -l
#PBS -N star_index
#PBS -l walltime=02:00:00
#PBS -l mem=50G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

## MANUAL SETUP BELOW

# >> SETUP: Specify Bowtie2 executable location
BT2DIR=/home/stewarz2/various_programs/bowtie2-2.4.5

# >> SETUP: Specify FASTA location
FASTADIR=/home/stewarz2/anemones/cassie/transcriptome/transcriptomes/trinity-denovo
FASTAFILE= #TBD

# >> SETUP: Specify computational resources
CPUS=1

## MANUAL SETUP END

## RUN PROGRAM
# STEP 1: Index the FASTA file where it resides
${BT2DIR}/bowtie2-build --threads ${CPUS} \
	${FASTADIR}/${FASTAFILE} \
	${FASTADIR}/${FASTAFILE}
