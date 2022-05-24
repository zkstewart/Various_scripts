#!/bin/bash -l
#PBS -N sam2bam
#PBS -l walltime=04:00:00
#PBS -l mem=20G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

## SETUP: Specify computational resources for qsub script
MEM=20

## AUTO SETUP: Derive the per-thread memory for samtools
BAMTOOLSMEM=$(echo "$(printf "%.0f\n" $(echo "(${MEM}*1000)"|bc -l))")

# STEP 1: Run samtools sort | bamtools index | samtools flagstat pipeline
for file in *.sam; do BASENAME=$(basename ${file}); PREFIX=${BASENAME%%.sam}; samtools view -b $file > ${PREFIX}.bam; bamtools sort -in ${PREFIX}.bam -out ${PREFIX}.sorted.bam -mem ${BAMTOOLSMEM}; bamtools index -in ${PREFIX}.sorted.bam; samtools flagstat ${PREFIX}.sorted.bam > ${PREFIX}.sorted.flagstat; done
