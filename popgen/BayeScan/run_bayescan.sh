#!/bin/bash -l

#PBS -N bayescan
#PBS -l ncpus=8
#PBS -l walltime=48:00:00
#PBS -l mem=80G

cd $PBS_O_WORKDIR

# START SETUP
# > START MANUAL SPECIFICATION
# >> Specify BayeScan exe location
BAYESCANDIR=/home/stewarz2/various_programs/BayeScan2.1/binaries
BAYESCAN=BayeScan2.1_linux64bits

# >> Specify computational parameters
CPUS=8

# >> Specify input file
INFILE=btrys06.final.geste

# > END MANUAL SPECIFICATION
# END SETUP

# START SCRIPT
# > STEP 1: Get our input files argument
${BAYESCANDIR}/${BAYESCAN} -threads ${CPUS} ${INFILE}

# END SCRIPT

