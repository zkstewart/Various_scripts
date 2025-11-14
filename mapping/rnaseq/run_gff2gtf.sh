#!/bin/bash -l
#PBS -N gff2gtf
#PBS -l walltime=01:00:00
#PBS -l mem=10G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

conda activate perl5
agat_convert_sp_gff2gtf.pl --gff manindi_flc.gff3 -o manindi_flc.gtf

