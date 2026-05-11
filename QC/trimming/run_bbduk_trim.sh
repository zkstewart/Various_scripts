#!/bin/bash -l
#PBS -N trim_reads
#PBS -l walltime=06:00:00
#PBS -l mem=35G
#PBS -l ncpus=4
#PBS -J 1-X

cd $PBS_O_WORKDIR

conda activate readqc # conda create -n readqc bioconda::bbmap

####

# Specify reads dir
READSDIR=/work/ePGL/sequencing/dna/illumina/citrus/COY18436_Dyneale/prepared_reads
R1SUFFIX=_1.fq.gz
R2SUFFIX=_2.fq.gz

#  Specify computational parameters
CPUS=4

####

# STEP 0: Locate the adapters.fa file
ADAPTERS=$(find ~/anaconda3/envs/readqc -type f -path "*/bbmap/resources/adapters.fa")

# STEP 1: Find file prefixes
declare -a RNAFILES
i=0
for f in ${READSDIR}/*${R1SUFFIX}; do
    RNAFILES[${i}]=$(echo "${f%%${R1SUFFIX}}");
    i=$((i+1));
done

# STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
FILEPREFIX=${RNAFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${FILEPREFIX})

# STEP 3: Run bbduk for trimming
if [[ ! -f ${BASEPREFIX}.ok ]]; then
    bbduk.sh -Xmx10g in1=${FILEPREFIX}${R1SUFFIX} in2=${FILEPREFIX}${R2SUFFIX} \
        out1=${BASEPREFIX}.trimmed_1P.fq.gz out2=${BASEPREFIX}.trimmed_2P.fq.gz \
        ref=${ADAPTERS} threads=${CPUS} \
        ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=5 minlength=25 tpe tbo && touch ${BASEPREFIX}.ok;
fi;
