#!/bin/bash -l
#PBS -N bwamem
#PBS -l walltime=12:00:00
#PBS -l mem=55G
#PBS -l ncpus=8
#PBS -J 1-X

cd $PBS_O_WORKDIR

####

# Specify reference genome FASTA location
GENOMEFASTA=/work/ePGL/genomes/mango/indica/CATAS_Mindica_2.1/alphonso_catas_2.1.fasta

# Specify reads dir
READSDIR=/work/ePGL/sequencing/rna/rnaseq/mango/NGS_526_Steph/trimmed_reads
R1SUFFIX=.trimmed_1P.fq.gz
R2SUFFIX=.trimmed_2P.fq.gz

# Specify computational parameters
CPUS=8

####

# STEP 1: Find read file prefixes
declare -a READFILES
i=0
for f in ${READSDIR}/*${R1SUFFIX}; do
    READFILES[${i}]=$(echo "${f%%${R1SUFFIX}}");
    i=$((i+1));
done

# STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
FILEPREFIX=${READFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${FILEPREFIX})

# STEP 3: Format readgroup
RG="'@RG\\tID:${BASEPREFIX}\\tSM:${BASEPREFIX}\\tPL:illumina\\tLB:lib1\\tPU:unit1'"

# STEP 4: Run bwa mem
cmd=$(echo "bwa mem -t ${CPUS} -R $(echo ${RG}) ${GENOMEFASTA} ${FILEPREFIX}${R1SUFFIX} ${FILEPREFIX}${R2SUFFIX} > ${BASEPREFIX}.sam")
eval $cmd
