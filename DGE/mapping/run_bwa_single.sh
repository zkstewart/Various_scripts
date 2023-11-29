#!/bin/bash -l
#PBS -N bwa
#PBS -l walltime=12:00:00
#PBS -l mem=60G
#PBS -l ncpus=4

cd $PBS_O_WORKDIR

####

# Specify BWA executable location
BWADIR=/home/stewarz2/various_programs/bwa

# Specify reads location
READSDIR=/home/stewarz2/banana_group/852_annotation/assembly/reads
R1=MAL-S_R1.fastq.gz
R2=MAL-S_R2.fastq.gz

# Specify genome file location
GENDIR=/home/stewarz2/banana_group/852_annotation/contigs
GENFILE=418_and_939.fasta

# Specify output prefix
PREFIX=MAL-S

# >> SETUP: Specify computational resources
CPUS=4

####

# STEP 1: Run BWA mapping
${BWADIR}/bwa mem -t ${CPUS} \
        ${GENDIR}/${GENFILE} \
        ${READSDIR}/${R1} ${READSDIR}/${R2} \
        > ${PREFIX}.sam

# STEP 2: Convert to BAM, sort, index, and flagstat for QC
samtools view --with-header -b ${PREFIX}.sam > ${PREFIX}.bam
samtools sort -@ 4 -O bam -o ${PREFIX}.sorted.bam -m 5G ${PREFIX}.bam
bamtools index -in ${PREFIX}.sorted.bam
samtools flagstat ${PREFIX}.sorted.bam > ${PREFIX}.sorted.flagstat
