#!/bin/bash -l
#PBS -N filtbwa
#PBS -l walltime=12:00:00
#PBS -l mem=55G
#PBS -l ncpus=8
#PBS -J 1-X

cd $PBS_O_WORKDIR

####

# Specify the HapHiC environment and install location
CONDAENV=haphic
HAPHICDIR=/home/stewarz2/various_programs/HapHiC

# Specify reads dir
READSDIR=/scratch/stewarz2/phase_variants/trimmed_hic
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

# STEP 3: Filter and convert to BAM
samtools view -h ${BASEPREFIX}.sam | samblaster | samtools view - -@ ${CPUS} -S -h -b -F 3340 -o ${BASEPREFIX}.hic.bam

# STEP 4: Run HapHiC-specific filtration
conda activate ${CONDAENV}
${HAPHICDIR}/utils/filter_bam ${BASEPREFIX}.hic.bam 1 --nm 3 --threads ${CPUS} | samtools view - -b -@ ${CPUS} -o ${BASEPREFIX}.haphic.bam
