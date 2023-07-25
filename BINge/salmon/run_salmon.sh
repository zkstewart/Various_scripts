#!/bin/bash -l
#PBS -N salmon
#PBS -l walltime=12:00:00
#PBS -l mem=60G
#PBS -l ncpus=8
#PBS -J 1-24

cd $PBS_O_WORKDIR

## MANUAL SETUP BELOW

# >> SETUP: Specify salmon executable location
SALMONDIR=/home/stewarz2/various_programs/salmon/salmon-1.9.0_linux_x86_64/bin

# >> SETUP: Specify salmon index location
INDEXDIR=/home/stewarz2/banana_group/qcav_dge/binge/salmon
INDEXFILE=qcav_dge_reference.fasta_salmon_index

# >> SETUP: Specify reads location & file suffix
## For the suffix, it's assumed that just prior to the given string there
## is the 1 / 2 suffix differentiating forward / reverse reads
READSDIR=/home/stewarz2/banana_group/qcav_dge/bbduk
R1SUFFIX=.trimmed_1P.fq
R2SUFFIX=.trimmed_2P.fq

# >> SETUP: Specify computational resources
CPUS=8

## MANUAL SETUP END

## RUN PROGRAM
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
    --index ${INDEXDIR}/${INDEXFILE} \
    --libType A \
    -1 ${PREFIX}${R1SUFFIX} \
    -2 ${PREFIX}${R2SUFFIX} \
    -o ${BASENAME}.out
