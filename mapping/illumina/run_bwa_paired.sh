#!/bin/bash -l
#PBS -N bwamem
#PBS -l walltime=12:00:00
#PBS -l mem=55G
#PBS -l ncpus=8
#PBS -J 1-30

cd $PBS_O_WORKDIR

####

# Specify reference genome
GENOMEDIR=/home/stewarz2/citrus/andrew_miles/powerpole/mapping/reference
GENOME=Power_both_polished.fasta

# Specify reads dir
READSDIR=/home/stewarz2/citrus/andrew_miles/powerpole/reads/dna/NGS_612
R1SUFFIX=.trimmed_1P.fq.gz
R2SUFFIX=.trimmed_2P.fq.gz

# Specify computational parameters
CPUS=8

####

# STEP 1: Find RNAseq file prefixes
declare -a RNAFILES
i=0
for f in ${READSDIR}/*${R1SUFFIX}; do
    RNAFILES[${i}]=$(echo "${f%%${R1SUFFIX}}");
    i=$((i+1));
done

# STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
FILEPREFIX=${RNAFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${FILEPREFIX})

# STEP 3: Format readgroup
RG="'@RG\\tID:${BASEPREFIX}\\tSM:${BASEPREFIX}\\tPL:illumina\\tLB:lib612\\tPU:unit1'"

# STEP 4: Run bwa mem
cmd=$(echo "bwa mem -t ${CPUS} -R $(echo ${RG}) ${GENOMEDIR}/${GENOME} ${FILEPREFIX}${R1SUFFIX} ${FILEPREFIX}${R2SUFFIX} > ${BASEPREFIX}.sam")
eval $cmd
