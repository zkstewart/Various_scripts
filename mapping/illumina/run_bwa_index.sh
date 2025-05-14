#!/bin/bash -l
#PBS -N bwa_index
#PBS -l walltime=01:30:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify BWA executable location
BWADIR=/home/stewarz2/various_programs/bwa

# Specify FASTA location
FASTADIR=/home/stewarz2/plant_group/ted/genome
FASTAFILE=CMJ.v1.0.gene.model_isos.trans

####

${BWADIR}/bwa index ${FASTADIR}/${FASTAFILE}
