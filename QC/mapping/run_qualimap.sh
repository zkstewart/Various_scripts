#!/bin/bash -l
#PBS -N qualimap
#PBS -l walltime=03:30:00
#PBS -l mem=15G
#PBS -l ncpus=1
#PBS -J 1-133

cd $PBS_O_WORKDIR

conda activate qualimap

####

BAMSUFFIX=.sorted.bam
MEM=10G

####

declare -i count
count=0
for file in *${BAMSUFFIX}; do
    BASENAME=$(basename ${file});
    PREFIX=${BASENAME%%${BAMSUFFIX}};
    count=($count+1);
    if [[ $count -eq ${PBS_ARRAY_INDEX} ]] ; then
        if [[ ! -f ${PREFIX}_bamqc/genome_results.txt  ]]; then
            qualimap bamqc -bam ${file} -outdir ${PREFIX}_bamqc --java-mem-size=${MEM};
        fi
    fi
done
