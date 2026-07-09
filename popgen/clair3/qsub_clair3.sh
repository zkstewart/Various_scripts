#!/bin/bash -l
#PBS -N clair3
#PBS -l walltime=48:00:00
#PBS -l mem=60G
#PBS -l ncpus=4
#PBS -J 1-15

cd $PBS_O_WORKDIR

conda activate clair3

####

# Specify the location of the genome FASTA
GENOME=/work/ePGL/genomes/citrus/murcott/henry_upuli_2025/murcott_hap1.fasta ## Make sure this is samtools faidx indexed first!

# Specify the location of the mapped BAM files
MAPDIR=/scratch/stewarz2/phase_variants/mapping/nanopore
BAMSUFFIX=.sorted.bam

# Specify the location of the model to be used
MODELDIR=/home/stewarz2/anaconda3/envs/clair3/bin/models
MODEL=r1041_e82_400bps_sup_v520

# Specify the computational resources to use
CPUS=4

#################################

# STEP 0: Get our file list
declare -a BAMFILES
i=0
for f in ${MAPDIR}/*${BAMSUFFIX}; do
    BAMFILES[${i}]=$(echo "${f}");
    i=$((i+1));
done

# STEP 1: Establish this subjob's variables
declare -i index
index=${PBS_ARRAY_INDEX}-1

INPUTFILE=${BAMFILES[${index}]}
PREFIX=$(basename ${INPUTFILE} ${BAMSUFFIX})

# STEP 2: Run Clair3
run_clair3.sh \
	--bam_fn=${INPUTFILE} \
	--ref_fn=${GENOME} \
	--model_path=${MODELDIR}/${MODEL} \
	--threads=${CPUS} \
	--platform="ont" \
	--output=${PREFIX} \
	--enable_phasing \
    --include_all_ctgs
