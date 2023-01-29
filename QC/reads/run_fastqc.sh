#!/bin/bash -l
#PBS -N fqc
#PBS -l walltime=24:00:00
#PBS -l mem=25G
#PBS -l ncpus=4
#PBS -J 1-12

cd $PBS_O_WORKDIR

### MANUAL SETUP BELOW
## SETUP: Load modules
module load fastqc/0.11.7-java-1.8.0_92

## SETUP: Specify RNAseq reads dir
READSDIR=/home/stewarz2/micro_group/carrie/reads
SUFFIX=.fastq.gz

## SETUP: Specify computational parameters
CPUS=4
# SETUP END
### MANUAL SETUP END


#####


# RUN START
## STEP 1: Find RNAseq file prefixes
declare -a RNAFILES
i=0
for f in ${READSDIR}/*1${SUFFIX}; do
    RNAFILES[${i}]=$(echo "${f%%1${SUFFIX}}");
    i=$((i+1));
done

## STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
FILEPREFIX=${RNAFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${FILEPREFIX} _) # strip any _ trailing characters

## STEP 3: Run FASTQC
mkdir -p ${BASEPREFIX}
fastqc -o ${BASEPREFIX} -t ${CPUS} ${FILEPREFIX}1${SUFFIX} ${FILEPREFIX}2${SUFFIX}
