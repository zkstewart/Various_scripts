#!/bin/bash -l
#PBS -N fb_thr
#PBS -l walltime=04:00:00
#PBS -l mem=5G
#PBS -l ncpus=12
#PBS -J 1-128

cd $PBS_O_WORKDIR

#################################

# Load the GCC that Freebayes was installed with (10.3.0)
module unload gcc/4.9.3-2.25
module unload gcccore/4.9.3
module unload binutils/2.25-gcccore-4.9.3
module unload zlib/1.2.8-foss-2016a
module load gcc/10.3.0
module unload libxml2/2.9.3-foss-2016a
module load libxml2/2.9.10-gcccore-10.3.0

# Specify things needed to run this
GENOMEDIR=/home/stewarz2/flies/chapa_2021/genome
GENOME=btrys06_freeze2.rename.fasta

MAPDIR=/home/stewarz2/flies/genome_based_2022/original/map/Parental_Selected

CPUS=12
NUMREGIONS=500 # number of chunks to parallelise to

#################################

# > STEP 1: Get our file list
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

# > STEP 5: Run freebayes with threading based on regions
freebayes-parallel ${PREFIX}.${NUMREGIONS}.regions ${CPUS} -f ${GENOMEDIR}/${GENOME} ${MAPDIR}/${INPUTFILE} > ${PREFIX}.vcf
