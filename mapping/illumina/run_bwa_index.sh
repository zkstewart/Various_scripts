#!/bin/bash -l
#PBS -N bwa_index
#PBS -l walltime=01:30:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify BWA executable location
BWADIR=/home/stewarz2/various_programs/bwa

# Specify reference genome FASTA location
GENOMEFASTA=/work/ePGL/genomes/mango/indica/CATAS_Mindica_2.1/alphonso_catas_2.1.fasta

####

${BWADIR}/bwa index ${GENOMEFASTA}
