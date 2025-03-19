#!/bin/bash -l
#PBS -N liftoff418to939
#PBS -l walltime=01:00:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

#################################

conda activate liftoff

# Specify the location of the input files
REFDIR=/home/stewarz2/banana_group/852_annotation/contigs
REFFILE=418.indelmodified.v2.fasta

GFFDIR=/home/stewarz2/banana_group/852_annotation/final_annotation
GFF=Ma418_models.v2.gff3

TARGETDIR=/home/stewarz2/banana_group/852_annotation/contigs
TARGETFILE=939.indelmodified.fasta

# Specify a prefix for output files
PREFIX=418_v2_to_939_v1

#################################

# Run liftoff
liftoff -g ${GFFDIR}/${GFF} \
        -o ${PREFIX}.gff3 \
        -copies -sc 0.9 \
        -polish \
        ${TARGETDIR}/${TARGETFILE} ${REFDIR}/${REFFILE}
