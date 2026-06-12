#!/bin/bash -l
#PBS -N tar719
#PBS -l walltime=12:00:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

READSDIR=downloaded
PREFIX=NGS719_Twist_test

####

tar --create --dereference --file=${PREFIX}.tar ${READSDIR}
md5sum ${PREFIX}.tar > ${PREFIX}.tar.md5
