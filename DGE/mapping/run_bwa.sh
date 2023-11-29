#!/bin/bash -l
#PBS -N bwa
#PBS -l walltime=12:00:00
#PBS -l mem=60G
#PBS -l ncpus=4
#PBS -J 1-12

cd $PBS_O_WORKDIR

####

# Specify BWA executable location
BWADIR=/home/stewarz2/various_programs/bwa

# Specify reads location
## For the suffix, it's assumed that just prior to the given string there
## is the 1 / 2 suffix differentiating forward / reverse reads
READSDIR=/home/stewarz2/banana_group/852_annotation/assembly/reads
SUFFIX=.fq.gz

# Specify genome file location
GENDIR=/home/stewarz2/banana_group/852_annotation/contigs
GENFILE=418_and_939.fasta

# Specify computational resources
CPUS=4

####

# STEP 1: Locate all read prefixes for mapping
declare -a RNAPREFIXES
i=0
for f in ${READSDIR}/*1${SUFFIX}; do
    RNAPREFIXES[${i}]=$(echo "${f%%1${SUFFIX}}");
    i=$((i+1));
done

# STEP 2: Handle batch submission variables
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
PREFIX=${RNAPREFIXES[${ARRAY_INDEX}]}
BASENAME=$(basename ${PREFIX} _) # strip any _ trailing characters

# STEP 3: Run BWA mapping
${BWADIR}/bwa mem -t ${CPUS} \
        ${GENDIR}/${GENFILE} \
        ${PREFIX}1${SUFFIX} ${PREFIX}2${SUFFIX} \
        > ${BASENAME}.sam

# STEP 4: Convert to BAM, sort, index, and flagstat for QC
samtools view --with-header -b ${BASENAME}.sam > ${BASENAME}.bam
samtools sort -@ 4 -O bam -o ${BASENAME}.sorted.bam -m 5G ${BASENAME}.bam
bamtools index -in ${BASENAME}.sorted.bam
samtools flagstat ${BASENAME}.sorted.bam > ${BASENAME}.sorted.flagstat
