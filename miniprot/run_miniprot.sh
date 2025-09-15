#!/bin/bash -l
#PBS -N miniprot
#PBS -l walltime=04:00:00
#PBS -l mem=25G
#PBS -l ncpus=12

cd $PBS_O_WORKDIR

####

# Specify the location of the BINge directory
MINIPROTDIR=/home/stewarz2/various_programs/miniprot

# Specify the locations of the genome and annotation files
GENOMEDIR=/home/stewarz2/plant_group/anuradha/upr_annotation/citrus_genomes
GENOME=australis_genome.fasta

PROTEINS=/home/stewarz2/plant_group/anuradha/upr_annotation/revised_models

# Specify computational parameters
THREADS=12

# Specify the output file prefix
OUTPREFIX=binge_results_update

####

# STEP 1: Index genome
${MINIPROTDIR}/miniprot -t${THREADS} -d ${GENOMEDIR}/${GENOME}.mpi ${GENOMEDIR}/${GENOME}

# STEP 2: Run alignment
${MINIPROTDIR}/miniprot -Iut${THREADS} --gff ${GENOMEDIR}/${GENOME}.mpi ${PROTEINS} > ${OUTPREFIX}.gff3
