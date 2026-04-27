#!/bin/bash -l
#PBS -N sam2bam
#PBS -l walltime=06:00:00
#PBS -l mem=35G
#PBS -l ncpus=1
#PBS -J 1-X

cd $PBS_O_WORKDIR

####

MEM=30

####

declare -i count
count=0
for file in *.sam; do
    BASENAME=$(basename ${file});
    PREFIX=${BASENAME%%.sam};
    count=($count+1);
    if [[ $count -eq ${PBS_ARRAY_INDEX} ]] ; then
        if [[ ! -f ${PREFIX}.bam  ]]; then
            samtools view -b ${file} > ${PREFIX}.bam;
        fi
        if [[ ! -f ${PREFIX}.sorted.bam ]]; then
            samtools sort -@ 1 -O bam -o ${PREFIX}.sorted.bam -m ${MEM}G ${PREFIX}.bam;
        fi
        if [[ ! -f ${PREFIX}.sorted.bam.bai ]]; then
            samtools index ${PREFIX}.sorted.bam;
        fi
        if [[ ! -f ${PREFIX}.sorted.flagstat ]]; then
            samtools flagstat ${PREFIX}.sorted.bam > ${PREFIX}.sorted.flagstat;
        fi
    fi
done
