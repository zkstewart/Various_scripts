#!/bin/bash -l
#PBS -N fbayes
#PBS -l walltime=00:10:00
#PBS -l mem=5G
#PBS -l ncpus=1
#PBS -J 1-159

cd $PBS_O_WORKDIR

#################################

# Specify the location of the freebayes executable
FBEXE=/home/stewarz2/various_programs/freebayes/freebayes-1.3.6-linux-amd64-static

# Specify the location of the genome FASTA
GENOMEDIR=/home/stewarz2/daniel/genome
GENOME=btrys06_chr1.fasta

# Specify the location of the mapped BAM files
MAPDIR=/home/stewarz2/daniel/rnaseq/chr1_map

# Specify the suffix that identifies mapped BAM files
SUFFIX=.sorted.md.bam

#################################

# > STEP 1: Get our file list
declare -a BAMFILES
i=0
for f in ${MAPDIR}/*${SUFFIX}; do
    BAMFILES[${i}]=$(echo "${f}");
    i=$((i+1));
done

# > STEP 2: Get our array index
declare -i index
index=${PBS_ARRAY_INDEX}-1

# > STEP 3: Get our file for analysis
INPUTFILE=${BAMFILES[${index}]}

# > STEP 4: Get our output file prefix
PREFIX=$(basename ${INPUTFILE} ${SUFFIX})

# > STEP 5: Run freebayes
${FBEXE} -f ${GENOMEDIR}/${GENOME} ${INPUTFILE} > ${PREFIX}.vcf

