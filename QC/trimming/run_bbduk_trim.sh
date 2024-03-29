#!/bin/bash -l
#PBS -N trim_bac
#PBS -l walltime=06:00:00
#PBS -l mem=35G
#PBS -l ncpus=4
#PBS -J 1-2

cd $PBS_O_WORKDIR

### MANUAL SETUP BELOW
## SETUP: Load modules
module load java/1.8.0_92

## SETUP: Specify bbmap location
BBDIR=/home/stewarz2/various_programs/bbmap

## SETUP: Specify RNAseq reads dir
READSDIR=/home/stewarz2/flies/mitch/prepared_reads
R1SUFFIX=_1.fq.gz
R2SUFFIX=_2.fq.gz

## SETUP: Specify computational parameters
CPUS=4
# SETUP END
### MANUAL SETUP END


#####


# RUN START
## STEP 1: Find RNAseq file prefixes
declare -a RNAFILES
i=0
for f in ${READSDIR}/*${R1SUFFIX}; do
    RNAFILES[${i}]=$(echo "${f%%${R1SUFFIX}}");
    i=$((i+1));
done

## STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
FILEPREFIX=${RNAFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${FILEPREFIX})

## STEP 3: Run bbduk for trimming
${BBDIR}/bbduk.sh -Xmx10g in1=${FILEPREFIX}${R1SUFFIX} in2=${FILEPREFIX}${R2SUFFIX} \
	out1=${BASEPREFIX}.trimmed_1P.fq out2=${BASEPREFIX}.trimmed_2P.fq \
	ref=${BBDIR}/resources/adapters.fa threads=${CPUS} \
	ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=5 minlength=25 tpe tbo
