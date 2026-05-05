#!/bin/bash -l
#PBS -N tarUMI
#PBS -l walltime=48:00:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

READSDIR=downloaded
PREFIX=UMI_citrus

####

tar --create --file=${PREFIX}.tar ${READSDIR}
md5sum ${PREFIX}.tar > ${PREFIX}.tar.md5
