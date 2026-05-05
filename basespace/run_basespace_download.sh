#!/bin/bash -l
#PBS -N fixedCOY18436
#PBS -l walltime=64:00:00
#PBS -l mem=35G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify basespace executable location
BSDIR=/home/stewarz2/various_programs/basespace

# Specify project ID
## Find this with: ${BSDIR}/bs projects list
PROJECTID=444155256

# Specify output directory
OUTDIR=downloaded

####

# Run basespace download
${BSDIR}/bs download project -i ${PROJECTID} -o ${OUTDIR}
