#!/bin/bash -l
#PBS -N k2tabulate
#PBS -l walltime=00:10:00
#PBS -l mem=10G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

conda activate kraken2 # should have bioconda::krakentools in this environment

####

# Specify the Various_scripts repository location
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the prefix for output file generation
PREFIX=NGS_526

####

combine_kreports.py -r */*_report.txt -o ${PREFIX}_kraken.tsv
python ${VARSCRIPTDIR}/kraken/kreport_reformat.py -i ${PREFIX}_kraken.tsv -o ${PREFIX}_kraken.reformat.tsv
