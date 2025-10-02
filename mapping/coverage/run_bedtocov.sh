#!/bin/bash -l
#PBS -N bed2cov
#PBS -l walltime=24:00:00
#PBS -l mem=40G
#PBS -l ncpus=1
#PBS -J 1-30

cd $PBS_O_WORKDIR

####

# Specify the genome lengths file location
GENOMEDIR=/home/stewarz2/citrus/andrew_miles/powerpole/mapping/reference
GENOMELENGTH=Power_both_polished.fasta.lengths.tsv

# Specify the suffix to identify files by
BEDSUFFIX=.bed

####

declare -i count
count=0
for file in *${BEDSUFFIX}; do
    BASENAME=$(basename ${file});
    PREFIX=${BASENAME%%${BEDSUFFIX}};
    count=($count+1);
    if [[ $count -eq ${PBS_ARRAY_INDEX} ]] ; then
        if [[ ! -f ${PREFIX}.cov  ]]; then
            bedtools genomecov -bg -i ${file} -g ${GENOMEDIR}/${GENOMELENGTH} > ${PREFIX}.cov;
        fi
    fi
done
