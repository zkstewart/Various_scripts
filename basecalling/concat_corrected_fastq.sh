#!/bin/bash -l
#PBS -N concat
#PBS -l walltime=03:00:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify the location of the Various_scripts and Genome_analysis_scripts folders
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts
GENSCRIPTDIR=/home/stewarz2/scripts/Genome_analysis_scripts

# Specify script behaviours
MINIMUM=5000
MULTILINE=60

####

# STEP 0: Set hard-coded variables
INDIR=step2
OUTDIR=output
mkdir -p ${OUTDIR}

# STEP 1: Concatenate individual FASTA files
if [[ ! -n "$(ls -A ${OUTDIR} 2>/dev/null)" ]]; then
    PREFIX=$(ls ${INDIR}/* | basename $(head -n 1) .min5000.fastq.block_0.fasta);
    python ${VARSCRIPTDIR}/fasta_concat.py \
        -i ${INDIR} \
        -o ${OUTDIR}/${PREFIX}_corrected_reads.fasta \
        --minimum ${MINIMUM} --multiline ${MULTILINE};
fi

# STEP 2: Obtain statistics for the reads
PREFIX=$(ls ${OUTDIR}/* | basename $(head -n 1) _corrected_reads.fasta)
python ${GENSCRIPTDIR}/genome_stats.py \
    -i ${OUTDIR}/${PREFIX}_corrected_reads.fasta \
    -o ${OUTDIR}/${PREFIX}_corrected_reads.stats
