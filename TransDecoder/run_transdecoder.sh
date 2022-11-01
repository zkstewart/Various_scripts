#!/bin/bash -l
#PBS -N trdecod
#PBS -l walltime=48:00:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# SETUP: Specify TransDecoder location
TRDECODEDIR=/home/stewarz2/various_programs/TransDecoder

# SETUP: Specify target file location
FASTADIR=/home/stewarz2/anemones/cassie/transcriptome/contaminant_removal/symbionts/downloads
FASTA=GSE72763_MI_min250_nr.fasta

####

# STEP 1: Run LongOrfs part of TransDecoder pipeline
${TRDECODEDIR}/TransDecoder.LongOrfs -t ${FASTADIR}/${FASTA}

# STEP 2: Run PRedict part of pipeline
${TRDECODEDIR}/TransDecoder.Predict -t ${FASTADIR}/${FASTA}

