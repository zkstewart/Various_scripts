#!/bin/bash -l
#PBS -N fb_prep
#PBS -l walltime=02:00:00
#PBS -l mem=120G
#PBS -l ncpus=1
#PBS -J 1-128

cd $PBS_O_WORKDIR

# Specify things needed to run this
GENOMEDIR=/home/stewarz2/flies/chapa_2021/genome
GENOME=btrys06_freeze2.rename.fasta

MAPDIR=/home/stewarz2/flies/genome_based_2022/original/map/Parental_Selected

NUMREGIONS=500 # number of chunks to parallelise to

# > STEP 1: Get our file list ## prefixes
cd ${MAPDIR}
declare -a BAMFILES
i=0
for f in *.sorted.bam; do
    BAMFILES[${i}]=$(echo "${f}");
    i=$((i+1));
done
cd $PBS_O_WORKDIR

# > STEP 2: Get our array index
declare -i index
index=${PBS_ARRAY_INDEX}-1

# > STEP 3: Get our file for analysis
INPUTFILE=${BAMFILES[${index}]}

# > STEP 4: Get our output file prefix
PREFIX=${INPUTFILE%%.sorted.bam}

# > STEP 5: Prep for freebayes in threaded mode by getting regions
bamtools coverage -in ${MAPDIR}/${INPUTFILE} | coverage_to_regions.py ${GENOMEDIR}/${GENOME}.fai ${NUMREGIONS} > ${PREFIX}.${NUMREGIONS}.regions
