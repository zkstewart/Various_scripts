#!/bin/bash -l
#PBS -N bam2bigwig
#PBS -l walltime=24:00:00
#PBS -l mem=40G
#PBS -l ncpus=1
#PBS -J 1-30

cd $PBS_O_WORKDIR

####

# Specify the bedGraphToBigWig exe location
BGTOBW_EXE=/home/stewarz2/various_programs/ucsc-twobit-toolkit/bedGraphToBigWig

# Specify the genome lengths file location
GENOMEDIR=/home/stewarz2/citrus/andrew_miles/powerpole/mapping/reference
GENOMELENGTH=Power_both_polished.fasta.lengths.tsv

# Specify the suffix to identify files by
BAMSUFFIX=.sorted.bam

####

declare -i count
count=0
for file in *${BAMSUFFIX}; do
    BASENAME=$(basename ${file});
    PREFIX=${BASENAME%%${BAMSUFFIX}};
    count=($count+1);
    if [[ $count -eq ${PBS_ARRAY_INDEX} ]] ; then
        if [[ ! -f ${PREFIX}.cov  ]]; then
            bedtools genomecov -bg -ibam ${file} -g ${GENOMEDIR}/${GENOMELENGTH} > ${PREFIX}.cov;
        fi
        if [[ ! -f ${PREFIX}.bigwig  ]]; then
            ${BGTOBW_EXE} ${PREFIX}.cov ${GENOMEDIR}/${GENOMELENGTH} ${PREFIX}.bigwig;
        fi
    fi
done
