#!/bin/bash -l
#PBS -N fptrim
#PBS -l walltime=01:30:00
#PBS -l mem=20G
#PBS -l ncpus=2
#PBS -J 1-10

cd $PBS_O_WORKDIR

module load GCC/13.2.0
module load fastp/0.23.4

####

# Specify RNAseq reads dir
READSDIR=/work/ePGL/sequencing/rna/rnaseq/macadamia/flowering_2021/prepared_reads
R1SUFFIX=_1.fq.gz
R2SUFFIX=_2.fq.gz

# Specify computational parameters
CPUS=2

#####

# STEP 1: Find RNAseq file prefixes
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

# STEP 3: Run fastp trimming
fastp --thread ${CPUS} \
    -i ${FILEPREFIX}${R1SUFFIX} -I ${FILEPREFIX}${R2SUFFIX} \
    -o ${BASEPREFIX}.trimmed_1P.fq.gz -O ${BASEPREFIX}.trimmed_2P.fq.gz \
    --detect_adapter_for_pe --trim_poly_x --cut_mean_quality 6 --length_required 25 \
    --html ${BASEPREFIX}.fastp.html --json ${BASEPREFIX}.fastp.json
