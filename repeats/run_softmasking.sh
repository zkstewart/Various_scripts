#!/bin/bash -l
#PBS -N softm
#PBS -l walltime=12:00:00
#PBS -l mem=30G
#PBS -l ncpus=12

cd $PBS_O_WORKDIR

####

# Specify program directories
REPMASKDIR=/home/stewarz2/various_programs/RepeatMasker

# Specify genome file location
GENOMEFILE=/work/ePGL/genomes/citrus/eg/something.hap1.fasta

# Specify repeats library file
REPLIB=/work/ePGL/databases/repeats/Citrus_repeatedSequences_library_v1.fasta

# Specify computational parameters
CPUS=12

# Specify output file prefix
PREFIX=something.hap1

####

# STEP 0: Put program into $PATH (since RepeatMasker expects this)
export PATH="${REPMASKDIR}:${PATH}"

# STEP 1: Run softmasking
mkdir -p ${PREFIX}_softmask

${REPMASKDIR}/RepeatMasker \
    -pa ${CPUS} \
    -lib ${REPLIB} \
    -dir ${PREFIX}_softmask \
    -e ncbi -s -nolow -no_is -norna -xsmall \
    ${GENOMEFILE}
