#!/bin/bash -l
#PBS -N orthof
#PBS -l walltime=32:00:00
#PBS -l mem=80G
#PBS -l ncpus=14

cd $PBS_O_WORKDIR

conda activate orthofinder3

####

# Specify input file location
FASTASDIR=/home/stewarz2/citrus/andrew_miles/powerpole/orthofinder/run

# Specify program resources
CPUS=14

# Specify output directory
OUTDIR=manual_v6_orthofinder

####

orthofinder -t ${CPUS} -a ${CPUS} -f ${FASTASDIR} -o ${OUTDIR}
