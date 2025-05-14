#!/bin/bash -l
#PBS -N insilicoNorm
#PBS -l walltime=48:00:00
#PBS -l mem=500G
#PBS -l ncpus=24

cd $PBS_O_WORKDIR

module load Java/17.0.6

####

# Specify Trinity location
TRINITYDIR=/home/stewarz2/various_programs/trinityrnaseq-v2.15.1

# Specify trimmed RNAseq reads dir
READSDIR=/work/ePGL/sequencing/rna/rnaseq/citrus/NGS_588_James/trimmed_reads
R1SUFFIX=.trimmed_1P.fq.gz
R2SUFFIX=.trimmed_2P.fq.gz

# Specify computational parameters
CPUS=24
MEM=450G

#####

# RUN START
# STEP 1: Concatenate left and right reads
cat ${READSDIR}/*${R1SUFFIX} > left.fq.gz
gunzip left.fq.gz
cat ${READSDIR}/*${R2SUFFIX} > right.fq.gz
gunzip right.fq.gz

## STEP 2: Run insilico normalisation
${TRINITYDIR}/util/insilico_read_normalization.pl --seqType fq \
    --JM ${MEM} --max_cov 30 --left left.fq --right right.fq \
    --pairs_together --PARALLEL_STATS --CPU ${CPUS} 2>&1 >> Trinity.log
