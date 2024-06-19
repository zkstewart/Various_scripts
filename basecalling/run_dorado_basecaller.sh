#!/bin/bash -l
#PBS -N dorado
#PBS -l walltime=24:00:00
#PBS -l mem=20G
#PBS -l ncpus=1
#PBS -l ngpus=1
#PBS -l gputype=A100

cd $PBS_O_WORKDIR

####

# Specify the location of the dorado bin and lib folders
DORADODIR=/home/stewarz2/various_programs/dorado/dorado-0.7.1-linux-x64

# Specify directory containing pod5 files for basecalling
POD5DIR=/work/ePGL/resequencing/nanopore/citrus_murcott_06_2024/Murcott_run2_a/pod5

# Specify model/model-complex to use
MODEL=sup

# Specify output prefix and directory to write basecalled results to
OUTDIR=/work/ePGL/resequencing/nanopore/citrus_murcott_06_2024/Murcott_run2_a/bam_sup
OUTPREFIX=pete_murcott_run2a

####

# STEP 0: Make sure LIB path is set and outdir exists
export LD_LIBRARY_PATH=${DORADODIR}/lib:$LD_LIBRARY_PATH
mkdir -p ${OUTDIR}

# STEP 0.5: Format model value for inclusion in output file name
MODELSUFFIX=$(echo ${MODEL} | tr , - | tr @ _)

# STEP 1: Run dorado
${DORADODIR}/bin/dorado basecaller ${MODEL} ${POD5DIR}/ > ${OUTDIR}/${OUTPREFIX}.${MODELSUFFIX}.bam
