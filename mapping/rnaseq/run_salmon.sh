#!/bin/bash -l
#PBS -N salmon
#PBS -l walltime=02:30:00
#PBS -l mem=15G
#PBS -l ncpus=4
#PBS -J 1-X

cd $PBS_O_WORKDIR

####

# Specify salmon executable location
SALMONDIR=/home/stewarz2/various_programs/salmon/salmon-1.9.0_linux_x86_64/bin

# Specify salmon index location
INDEXDIR=/work/ePGL/genomes/citrus/murcott/henry_upuli_2025/murcott_hap1.exon.fasta_salmon_index

# Specify reads location & file suffix
## For the suffix, it's assumed that just prior to the given string there
## is the 1 / 2 suffix differentiating forward / reverse reads
READSDIR=/scratch/stewarz2/citrus/andrew_miles/dge_experiments/ngs713/trimmed_reads
R1SUFFIX=.trimmed_1P.fq
R2SUFFIX=.trimmed_2P.fq

# Specify computational resources
CPUS=4

####

# STEP 1: Locate all read prefixes for mapping
declare -a RNAPREFIXES
i=0
for f in ${READSDIR}/*${R1SUFFIX}; do
    RNAPREFIXES[${i}]=$(echo "${f%%${R1SUFFIX}}");
    i=$((i+1));
done

# STEP 2: Handle batch submission variables
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
PREFIX=${RNAPREFIXES[${ARRAY_INDEX}]}
BASENAME=$(basename ${PREFIX})

# STEP 3: Run mapping procedure for each sample
${SALMONDIR}/salmon quant --threads ${CPUS} \
                          --index ${INDEXDIR} \
                          --libType A \
                          -1 ${PREFIX}${R1SUFFIX} \
                          -2 ${PREFIX}${R2SUFFIX} \
                          -o ${BASENAME}.out
