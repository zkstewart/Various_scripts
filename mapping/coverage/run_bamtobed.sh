#!/bin/bash -l
#PBS -N bam2bed
#PBS -l walltime=24:00:00
#PBS -l mem=40G
#PBS -l ncpus=1
#PBS -J 1-30

cd $PBS_O_WORKDIR

####

BAMSUFFIX=.sorted.bam

####

declare -i count
count=0
for file in *${BAMSUFFIX}; do
    BASENAME=$(basename ${file});
    PREFIX=${BASENAME%%${BAMSUFFIX}};
    count=($count+1);
    if [[ $count -eq ${PBS_ARRAY_INDEX} ]] ; then
        if [[ ! -f ${PREFIX}.bed  ]]; then
            bedtools bamtobed -i ${file} > ${PREFIX}.bed;
        fi
        if [[ ! -f ${PREFIX}.bedpe  ]]; then
            bedtools bamtobed -bedpe -i ${file} > ${PREFIX}.bedpe;
        fi
    fi
done
