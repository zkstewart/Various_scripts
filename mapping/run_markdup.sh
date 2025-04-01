#!/bin/bash -l
#PBS -N markdup
#PBS -l walltime=06:00:00
#PBS -l mem=20G
#PBS -l ncpus=1
#PBS -J 1-X

cd $PBS_O_WORKDIR

module load picard/3.0.0-Java-17

####

declare -i count
count=0
for file in *.sorted.bam; do
    BASENAME=$(basename ${file});
    PREFIX=${BASENAME%%.sorted.bam};
    count=($count+1);
    if [[ $count -eq ${PBS_ARRAY_INDEX} ]] ; then
        if [[ ! -f ${PREFIX}.sorted.md.bam  ]]; then
            java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$file O=${PREFIX}.sorted.md.bam M=${PREFIX}.md_metrics.txt;
        fi
    fi
done
