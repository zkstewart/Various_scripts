#!/bin/bash -l
#PBS -N trim_hic
#PBS -l walltime=01:00:00
#PBS -l mem=35G
#PBS -l ncpus=4
#PBS -J 1-2

cd $PBS_O_WORKDIR

####

# Specify fastp location
FPDIR=/home/stewarz2/various_programs/fastp

# Specify reads dir
READSDIR=/work/ePGL/sequencing/dna/hic/citrus/NGS_651/prepared_reads
R1SUFFIX=_1.fq.gz
R2SUFFIX=_2.fq.gz

# Specify computational parameters
CPUS=4

####

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

# STEP 3: Run fastp for trimming
${FPDIR}/fastp -i ${FILEPREFIX}${R1SUFFIX} -I ${FILEPREFIX}${R2SUFFIX} \
    -o ${BASEPREFIX}.trimmed_1P.fq.gz -O ${BASEPREFIX}.trimmed_2P.fq.gz \
    --thread ${CPUS} --trim_poly_g
