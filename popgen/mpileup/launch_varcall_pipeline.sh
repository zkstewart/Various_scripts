#!/bin/bash -l

####

# Specify the location of the Various_scripts repository
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the location of the genome FASTA
GENOME=/home/stewarz2/plant_group/soh_kimura/genomes/clementina.genome.fasta

# Specify the location of the mapped BAM files
BAMDIR=/home/stewarz2/plant_group/soh_kimura/mapping/results

# Specify the suffix that identifies mapped BAM files
SUFFIX=.sorted.md.bam

# Specify the prefix for output files
PREFIX=UQCraigSohKimura

# Specify computational parameters
WALLTIME=48:00:00
MEM=60G

####

python ${VARSCRIPTDIR}/popgen/mpileup/variant_calling_pipeline.py \
    -f ${GENOME} -b ${BAMDIR} --bamSuffix ${SUFFIX} \
    --walltime ${WALLTIME} --mem ${MEM} --jobPrefix ${PREFIX} \
    --indels-cns
