#!/bin/bash -l
#PBS -N salmon_index
#PBS -l walltime=24:00:00
#PBS -l mem=65G
#PBS -l ncpus=8

cd $PBS_O_WORKDIR

####

# Specify salmon executable location
SALMONDIR=/home/stewarz2/various_programs/salmon/salmon-1.9.0_linux_x86_64/bin

# Specify FASTA location
FASTAFILE=murcott_hap1.exon.fasta

# Specify computational resources
CPUS=8

####

${SALMONDIR}/salmon index --threads ${CPUS} \
                          --transcripts ${FASTAFILE} \
                          --index ${FASTADIR}/${FASTAFILE}_salmon_index
