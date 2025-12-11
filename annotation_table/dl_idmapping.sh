#!/bin/bash -l
#PBS -N dl_mapping
#PBS -l walltime=48:00:00
#PBS -l mem=15G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

wget --no-verbose https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
gunzip idmapping_selected.tab.gz
