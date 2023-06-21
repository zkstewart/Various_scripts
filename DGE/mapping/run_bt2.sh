#!/bin/bash -l
#PBS -N bowtie2
#PBS -l walltime=04:00:00
#PBS -l mem=50G
#PBS -l ncpus=4
#PBS -J 1-12

cd $PBS_O_WORKDIR

## MANUAL SETUP BELOW

# >> SETUP: Specify Bowtie2 executable location
BT2DIR=/home/stewarz2/various_programs/bowtie2-2.4.5

# >> SETUP: Specify FASTA location
FASTADIR=/home/stewarz2/anemones/cassie/transcriptome/transcriptomes/trinity-denovo
FASTAFILE= #TBD

# >> SETUP: Specify reads location & file suffix
## For the suffix, it's assumed that just prior to the given string there
## is the 1 / 2 suffix differentiating forward / reverse reads
READSDIR=/home/stewarz2/anemones/cassie/transcriptome/trimmomatic
SUFFIX=P.trimmed.fq

# >> SETUP: Specify computational resources
CPUS=4

## MANUAL SETUP END

## RUN PROGRAM
# STEP 1: Locate all read prefixes for mapping
declare -a RNAPREFIXES
i=0
for f in ${READSDIR}/*1${SUFFIX}; do
    RNAPREFIXES[${i}]=$(echo "${f%%1${SUFFIX}}");
    i=$((i+1));
done

# STEP 2: Handle batch submission variables
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
PREFIX=${RNAPREFIXES[${ARRAY_INDEX}]}
BASENAME=$(basename ${PREFIX} _) # strip any _ trailing characters

# STEP 3: Setup directory for output files
mkdir -p ${BASENAME}
cd ${BASENAME}

# STEP 4: Run mapping procedure for each sample
${BT2DIR}/bowtie2 --threads ${CPUS} \
	--sensitive \
	-x ${FASTADIR}/${FASTAFILE} \
	-1 ${PREFIX}1${SUFFIX} \
	-2 ${PREFIX}2${SUFFIX} \
	-S ${BASENAME}.sam
