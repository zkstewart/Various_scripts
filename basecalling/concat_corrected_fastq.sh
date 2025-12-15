#!/bin/bash -l
#PBS -N concat
#PBS -l walltime=01:00:00
#PBS -l mem=30G
#PBS -l ncpus=1
#PBS -W depend=afterok:

cd $PBS_O_WORKDIR

####

# Specify the location of the Various_scripts folder
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify file prefix
PREFIX=13Q030

# Specify script behaviours
MINIMUM=5000
MULTILINE=60

####

# STEP 0: Set hard-coded variables
INDIR=step2
OUTDIR=output

# STEP 1: Run samtools to convert BAM to FASTQ
mkdir -p ${OUTDIR}
if [[ ! -f ${OUTDIR}/${PREFIX}_corrected_reads.fasta ]]; then
    python ${VARSCRIPTDIR}/fasta_concat.py -i ${INDIR} -o ${OUTDIR}/${PREFIX}_corrected_reads.fasta \
        --minimum ${MINIMUM} --multiline ${MULTILINE};
fi
