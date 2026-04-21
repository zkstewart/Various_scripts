#!/bin/bash -l
#PBS -N concat
#PBS -l walltime=03:00:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify the location of the Various_scripts folder
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify script behaviours
MINIMUM=5000
MULTILINE=60

####

# STEP 0: Set hard-coded variables
INDIR=step2
OUTDIR=output
mkdir -p ${OUTDIR}

# STEP 1: Derive the file prefix
EGFILE=$(ls ${INDIR}/* | basename $(head -n 1))
SAMPLENUM=$(echo ${EGFILE} | cut -d "_" -f3 | cut -d "." -f1)
PREFIX=NGS_662_${SAMPLENUM}

# STEP 2: Concatenate individual FASTA files
python ${VARSCRIPTDIR}/fasta_concat.py \
    -i ${INDIR} \
    -o ${OUTDIR}/${PREFIX}_corrected_reads.fasta \
    --minimum ${MINIMUM} --multiline ${MULTILINE};
