#!/bin/bash -l
#PBS -N fs_choosek
#PBS -l walltime=00:20:00
#PBS -l mem=5G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

# >> Load in relevant modules
conda activate py2
module load gsl/2.1-foss-2016a
export CFLAGS="-I/pkg/suse12/software/gsl/2.1-foss-2016a/include"
export LDFLAGS="-L/pkg/suse12/software/gsl/2.1-foss-2016a/lib"

# >> Specify location of fastStructure install
FSDIR=/home/stewarz2/old/various_programs/fastStructure

# >> Specify location of multiple K runs
RUNSDIR=/home/stewarz2/flies/chapa_2022/fastSTRUCTURE/analyse/simple

# >> Specify file prefix
PREFIX=btrys06_simple

# STEP 1: Run script
python ${FSDIR}/chooseK.py --input=${RUNSDIR}/${PREFIX}
