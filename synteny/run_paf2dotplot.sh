#!/bin/bash -l
#PBS -N paf2dotplot
#PBS -l walltime=03:30:00
#PBS -l mem=30G
#PBS -l ncpus=8

cd $PBS_O_WORKDIR

conda activate synteny

#################################

# Specify the location of the paf2dotplot file
P2DPDIR=/home/stewarz2/various_programs/paf2dotplot

# Specify the location of the input files
REF=/home/stewarz2/path/to/ref.fasta
QUERY=/home/stewarz2/path/to/query.fasta

# Specify how many CPUs to use
CPUS=8

# Specify prefix for outputs
PREFIX=outputPrefix

# Specify behavioural params
INCHES=15
MINREFLEN=1000000
MINQUERYLEN=400000
MINALNLEN=10000

#################################

# STEP 1: Align with minimap2
minimap2 -cx asm20 -t ${CPUS} \
    ${REF} ${QUERY} > ${PREFIX}.paf

# STEP 2: Use syri to find structural annotations between genomes
${P2DPDIR}/paf2dotplot.r -f -b \
    --plot-size=${INCHES} \
    --min-ref-len=${MINREFLEN} \
    --min-query-length=${MINQUERYLEN} \
    --min-alignment-length=${MINALNLEN} \
    -o ${PREFIX} \
     ${PREFIX}.paf
