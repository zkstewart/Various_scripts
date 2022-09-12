#!/bin/bash -l
#PBS -N orthoF_citrus
#PBS -l walltime=48:00:00
#PBS -l mem=40G
#PBS -l ncpus=8

cd $PBS_O_WORKDIR

module load mafft/7.305-foss-2016a-with-extensions

# Setup: Specify the location of OrthoFinder's python file
ORTHODIR=/home/stewarz2/various_programs/OrthoFinder

# Setup: Specify input file locations
FASTASDIR=/home/stewarz2/citrus/orthologs/orthofinder_fastas

# Setup: Manual specification of program resources
CPUS=8


# RUN PROGRAM
python ${ORTHODIR}/orthofinder.py -t ${CPUS} -a ${CPUS} -f ${FASTASDIR}
