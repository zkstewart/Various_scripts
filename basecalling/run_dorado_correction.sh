#!/bin/bash -l

JOBPREFIX=18Q053
MINREADLEN=5000

####

# Specify the location of the dorado bin and lib folders
DORADODIR=/home/stewarz2/various_programs/dorado/dorado-1.3.0-linux-x64

# Specify the location of the Various_scripts folder
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts 

# Specify FASTQ file for correction
FQDIR=fastq
FQFILE=${JOBPREFIX}.min${MINREADLEN}.fastq

# Specify output directory
OUTDIR=fastq_corrected

####

NUMBLOCKS=$(cat ${FQDIR}/numblocks.txt)

python ${VARSCRIPTDIR}/basecalling/dorado_correct_pipeline.py -f ${FQDIR}/${FQFILE} \
    -n ${NUMBLOCKS} -o ${OUTDIR} \
    --dorado ${DORADODIR}/bin/dorado \
    --jobPrefix ${JOBPREFIX} --memCPU 650G --walltimeCPU 48:00:00 --cpuCPU 22 --walltimeGPU 16:00:00 --memGPU 60G
