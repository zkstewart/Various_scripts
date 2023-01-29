#!/bin/bash -l
#PBS -N star
#PBS -l walltime=04:00:00
#PBS -l mem=100G
#PBS -l ncpus=4
#PBS -J 1-12

cd $PBS_O_WORKDIR

## MANUAL SETUP BELOW

# >> SETUP: Specify STAR executable location
STARDIR=/home/stewarz2/various_programs/STAR-2.7.10a/bin/Linux_x86_64_static

# >> SETUP: Specify reads location & file suffix
## For the suffix, it's assumed that just prior to the given string there
## is the 1 / 2 suffix differentiating forward / reverse reads
READSDIR=/home/stewarz2/plant_group/leena/analysis_2/rnaseq_reads
SUFFIX=.fq.gz

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
## We set up cmds in an array so we can conditionally add something to it
cmd=(${STARDIR}/STAR --runThreadN ${CPUS} \
        --genomeDir ${PBS_O_WORKDIR} \
        --readFilesIn ${PREFIX}1${SUFFIX} ${PREFIX}2${SUFFIX} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --quantMode TranscriptomeSAM GeneCounts)
if [[ ${SUFFIX} == *.gz ]] ; then
	cmd+=(--readFilesCommand zcat);
fi;
"${cmd[@]}"
