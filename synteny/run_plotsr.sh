#!/bin/bash -l
#PBS -N plotsrPipe
#PBS -l walltime=03:30:00
#PBS -l mem=30G
#PBS -l ncpus=8

cd $PBS_O_WORKDIR

conda activate synteny

#################################

## This script DOES require you to change stuff outside of the '####' block here.
## Specifically, step 3 will require you to specify the genomes.txt file appropriately.^S

# Specify the location of the input files
GEN1=/home/stewarz2/plant_group/juel/genome/glauca.fasta
GEN2=/home/stewarz2/plant_group/juel/genome/hindsii.fasta

# Specify how many CPUs to use
CPUS=8

# Specify prefix for outputs
PREFIX=glauca_comparison

#################################

# STEP 1: Align with minimap2
minimap2 -ax asm5 -t ${CPUS} --eqx ${GEN1} ${GEN2} \
 | samtools sort -O BAM - > ${PREFIX}.bam
samtools index ${PREFIX}.bam

# STEP 2: Use syri to find structural annotations between genomes
syri -c ${PREFIX}.bam -r ${GEN1} -q ${GEN2} -F B --prefix ${PREFIX}

# STEP 3: Generate genomes.txt file for plotsr to interpret
echo "#file	name	tags" > genomes.txt
echo "${GEN1}	glauca	lw:1.5" >> genomes.txt
echo "${GEN2}	hindsii	lw:1.5" >> genomes.txt

# STEP 4: Run plotsr
plotsr \
    --sr ${PREFIX}syri.out \
    --genomes genomes.txt \
    -o ${PREFIX}_plot.png
