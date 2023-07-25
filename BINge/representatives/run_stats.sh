#!/bin/bash -l
#PBS -N BINge_stats
#PBS -l walltime=01:00:00
#PBS -l mem=15G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

# Specify the location of the BINge directory
GENSCRIPTDIR=/home/stewarz2/scripts/Genome_analysis_scripts
PREFIX=qcav_representatives.BINge.filtered

python ${GENSCRIPTDIR}/genome_stats.py -i ${PREFIX}.aa -o ${PREFIX}.aa.stats
python ${GENSCRIPTDIR}/genome_stats.py -i ${PREFIX}.cds -o ${PREFIX}.cds.stats
python ${GENSCRIPTDIR}/genome_stats.py -i ${PREFIX}.fasta -o ${PREFIX}.fasta.stats
