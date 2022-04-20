#!/bin/bash -l
#PBS -N fs_simple
#PBS -l walltime=12:00:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -J 1-30

cd $PBS_O_WORKDIR

# >> Load in relevant modules
conda activate py2
module load gsl/2.1-foss-2016a
export CFLAGS="-I/pkg/suse12/software/gsl/2.1-foss-2016a/include"
export LDFLAGS="-L/pkg/suse12/software/gsl/2.1-foss-2016a/lib"

# >> Specify location of fastStructure install
FSDIR=/home/stewarz2/old/various_programs/fastStructure

# >> Specify input file location
INDIR=/home/stewarz2/flies/chapa_2022/format_conversions
INFILE=btrys06.final.split

# >> Specify output prefix
OUTPREF=btrys06


# STEP 1: Run structure.py with various K values
python ${FSDIR}/structure.py -K ${PBS_ARRAY_INDEX} --input=${INDIR}/${INFILE} --output=${OUTPREF}_simple --full --seed=100 --format=bed
