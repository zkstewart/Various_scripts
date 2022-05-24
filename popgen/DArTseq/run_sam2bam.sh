#!/bin/bash -l
#PBS -N sam2bam
#PBS -l walltime=02:00:00
#PBS -l mem=50G
#PBS -l ncpus=4

cd $PBS_O_WORKDIR

## SETUP: Specify computational resources for qsub script
CPUS=4
MEM=50

## AUTO SETUP: Derive the per-thread memory for samtools
SAMTOOLSTHREADMEM=$(echo "$(printf "%.0f\n" $(echo "(${MEM}*0.50)/${CPUS}"|bc -l))")

# STEP 1: Run samtools sort | bamtools index | samtools flagstat pipeline
for file in *.sam; do BASENAME=$(basename ${file}); PREFIX=${BASENAME%%.sam}; samtools sort -m ${SAMTOOLSTHREADMEM}G -@ ${CPUS} $file > ${PREFIX}.sorted.bam; bamtools index -in ${PREFIX}.sorted.bam; samtools flagstat ${PREFIX}.sorted.bam > ${PREFIX}.sorted.flagstat; done

