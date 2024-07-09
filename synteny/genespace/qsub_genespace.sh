#!/bin/bash -l
#PBS -N allGENESPACE
#PBS -l walltime=120:00:00
#PBS -l mem=100G
#PBS -l ncpus=12

cd $PBS_O_WORKDIR

conda activate synteny

####

Rscript run_genespace.R --no-save
