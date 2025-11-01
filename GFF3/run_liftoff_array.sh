#!/bin/bash -l
#PBS -N liftoff
#PBS -l walltime=06:00:00
#PBS -l mem=40G
#PBS -l ncpus=8
#PBS -J 1-26

cd $PBS_O_WORKDIR

conda activate liftoff

####

# Specify the location of the reference files
REFGENOME=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/references/concatenated/Busk.ref.bothhap.fasta
REFGFF3=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/references/concatenated/Busk.ref.bothhap.gff3

# Specify the directory of target files
TARGETDIR=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/hap1
TARGETSUFFIX=.fasta

# Specify computational parameters
CPUS=8

####

# STEP 1: Locate all target genome files
declare -a TARGETFILES
i=0
for f in ${TARGETDIR}/*${TARGETSUFFIX}; do
    TARGETFILES[${i}]=$(echo "${f%%${TARGETSUFFIX}}");
    i=$((i+1));
done

# STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
TARGETPREFIX=${TARGETFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${TARGETPREFIX})

# STEP 3: Run liftoff on target genome
liftoff -p ${CPUS} \
    -o ${BASEPREFIX}.gff3 -dir tmp_${BASEPREFIX} \
    -copies -sc 0.9 -polish \
    -g ${REFGFF3} ${TARGETPREFIX}${TARGETSUFFIX} ${REFGENOME};
