#!/bin/bash -l

JOBPREFIX=roberts_sample5
NUMBLOCKS=7

####

# Specify the location of the dorado bin and lib folders
DORADODIR=/home/stewarz2/various_programs/dorado/dorado-0.9.0-linux-x64

# Specify the location of the Various_scripts folder
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts 

# Specify FASTQ file for correction
FQDIR=/work/ePGL/resequencing/nanopore/glauca_shearing_tests/Pete_Roberts_sample5__08_10_2024/fastq
FQFILE=glauca_shearing_roberts_sample5.min5000.fastq

# Specify output directory
OUTDIR=fastq_corrected

####

python ${VARSCRIPTDIR}/basecalling/dorado_correct_pipeline.py -f ${FQDIR}/${FQFILE} \
    -n ${NUMBLOCKS} -o ${OUTDIR} \
    --dorado ${DORADODIR}/bin/dorado \
    --jobPrefix ${JOBPREFIX} --memCPU 650G --walltimeCPU 48:00:00 --cpuCPU 22 --walltimeGPU 16:00:00 --memGPU 60G
