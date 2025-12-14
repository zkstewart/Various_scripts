#!/bin/bash -l
#PBS -N subcorr
#PBS -l walltime=00:10:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -W depend=afterok:

cd $PBS_O_WORKDIR

####

# Specify parameters for this sample / basecalling run
PREFIX=18Q053
MINREADLEN=5000

# Specify the location of the dorado bin and lib folders
DORADODIR=/home/stewarz2/various_programs/dorado/dorado-1.3.0-linux-x64

# Specify the location of the Various_scripts folder
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts 

# Specify input and output read file locations
FQDIR=fastq
OUTDIR=fastq_corrected

####

# STEP 0: Derive file locations (following on from run_dorado_basecaller.sh)
FQFILE=${PREFIX}.min${MINREADLEN}.fastq
NUMBLOCKS=$(cat ${FQDIR}/numblocks.txt)

# STEP 1: Run the pipeliner script
python ${VARSCRIPTDIR}/basecalling/dorado_correct_pipeline.py -f ${FQDIR}/${FQFILE} \
    -n ${NUMBLOCKS} -o ${OUTDIR} \
    --dorado ${DORADODIR}/bin/dorado \
    --jobPrefix ${PREFIX} --memCPU 650G --walltimeCPU 48:00:00 --cpuCPU 22 --walltimeGPU 16:00:00 --memGPU 60G
