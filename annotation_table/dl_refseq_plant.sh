#!/bin/bash -l
#PBS -N dl_refseq
#PBS -l walltime=24:00:00
#PBS -l mem=20G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

wget --no-verbose https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER
wget --no-verbose ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/plant.*.protein.faa.gz

cat plant.*.protein.faa.gz > refseq_plant.protein.faa.gz
gunzip refseq_plant.protein.faa.gz
