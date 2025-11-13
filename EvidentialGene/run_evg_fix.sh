#!/bin/bash -l
#PBS -N evgfix
#PBS -l walltime=00:30:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify Various_scripts location
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify EvidentialGene run dir
#EVGRUNDIR=/home/stewarz2/citrus/andrew_miles/powerpole/transcriptome/toxin2/transcriptome/transcriptomes/evidentialgene/evgrun
EVGRUNDIR=../evgrun

####

python ${VARSCRIPTDIR}/EvidentialGene/fix_evg_missing_orfs.py \
    -drop ${EVGRUNDIR}/dropset \
    -fasta okay-okalt.fasta \
    -aa okay-okalt.aa \
    -cds okay-okalt.cds \
    -o okay-okalt.final
