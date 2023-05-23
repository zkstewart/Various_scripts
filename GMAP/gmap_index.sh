#!/bin/bash -l
#PBS -N gmapIndex
#PBS -l walltime=04:00:00
#PBS -l mem=20G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

# Specify the location of the genome file
GENDIR=/home/stewarz2/telmatactis
GENFILE=telmatactis_HGAP.arr4.pil2.fasta

#####

# > STEP 1: Format GMAP database
gmap_build -D ${GENDIR} -d ${GENFILE}.gmap ${GENDIR}/${GENFILE}
