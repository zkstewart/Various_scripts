#!/bin/bash -l
#PBS -N trim_reads
#PBS -l walltime=06:00:00
#PBS -l mem=35G
#PBS -l ncpus=4
#PBS -J 1-2

cd $PBS_O_WORKDIR

module load Java/17.0.6

####

# Specify bbmap location
BBDIR=/home/stewarz2/various_programs/bbmap

# Specify reads dir
READSDIR=/home/stewarz2/flies/mitch/prepared_reads
R1SUFFIX=_1.fq.gz
R2SUFFIX=_2.fq.gz

#  Specify computational parameters
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

# STEP 3: Run bbduk for trimming
${BBDIR}/bbduk.sh -Xmx10g in1=${FILEPREFIX}${R1SUFFIX} in2=${FILEPREFIX}${R2SUFFIX} \
	out1=${BASEPREFIX}.trimmed_1P.fq.gz out2=${BASEPREFIX}.trimmed_2P.fq.gz \
	ref=${BBDIR}/resources/adapters.fa threads=${CPUS} \
	ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=5 minlength=25 tpe tbo
