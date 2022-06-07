#!/bin/bash -l
#PBS -N markdup
#PBS -l walltime=04:00:00
#PBS -l mem=20G
#PBS -l ncpus=1
#PBS -J 1-73

cd $PBS_O_WORKDIR

module load atg/picard/2.2.2

# STEP 1: Run Picard MarkDuplicates
declare -i count
count=0
for file in *.sorted.bam; do
    BASENAME=$(basename ${file});
    PREFIX=${BASENAME%%.sorted.bam};
    count=($count+1);
    if [[ $count -eq ${PBS_ARRAY_INDEX} ]] ; then
        if [[ ! -f ${PREFIX}.sorted.md.bam  ]]; then
            picard MarkDuplicates I=$file O=${PREFIX}.sorted.md.bam M=${PREFIX}.md_metrics.txt;
        fi
    fi
done

