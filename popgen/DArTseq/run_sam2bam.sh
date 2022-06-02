#!/bin/bash -l
#PBS -N sam2bam
#PBS -l walltime=08:00:00
#PBS -l mem=20G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

## SETUP: Specify computational resources for qsub script
MEM=20

## AUTO SETUP: Derive the per-thread memory for samtools
BAMTOOLSMEM=$(echo "$(printf "%.0f\n" $(echo "(${MEM}*1000)"|bc -l))")

# STEP 1: Run samtools view | bamtools sort | bamtools index | samtools flagstat pipeline
for file in *.sam; do
    BASENAME=$(basename ${file});
    PREFIX=${BASENAME%%.sam};
    if [[ ! -f ${PREFIX}.bam  ]]; then
        samtools view -b $file > ${PREFIX}.bam;
    fi
    if [[ ! -f ${PREFIX}.sorted.bam ]]; then
        bamtools sort -in ${PREFIX}.bam -out ${PREFIX}.sorted.bam -mem ${BAMTOOLSMEM};
    fi
    if [[ ! -f ${PREFIX}.sorted.bam.bai ]]; then
        bamtools index -in ${PREFIX}.sorted.bam;
    fi
    if [[ ! -f ${PREFIX}.sorted.flagstat ]]; then
        samtools flagstat ${PREFIX}.sorted.bam > ${PREFIX}.sorted.flagstat;
    fi
done
