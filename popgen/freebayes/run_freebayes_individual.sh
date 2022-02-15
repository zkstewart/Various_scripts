#!/bin/bash -l
#PBS -N freebayes
#PBS -l walltime=12:00:00
#PBS -l mem=20G
#PBS -l ncpus=1
#PBS -J 1-186

cd $PBS_O_WORKDIR

# Specify things needed to run this
GENOMEDIR=/home/stewarz2/flies/chapa_2021/genome
GENOME=btrys06_freeze2.rename.fasta

MAPDIR=/home/stewarz2/flies/chapa_2022/map


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

# > STEP 4: Run freebayes
freebayes -f ${GENOMEDIR}/${GENOME} ${MAPDIR}/${INPUTFILE} > ${PREFIX}.vcf
