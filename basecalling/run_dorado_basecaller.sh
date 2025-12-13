#!/bin/bash -l
#PBS -N dorado
#PBS -l walltime=48:00:00
#PBS -l mem=60G
#PBS -l ncpus=1
#PBS -l ngpus=1
#PBS -l gpu_id=A100

cd $PBS_O_WORKDIR

module load zlib/1.3.1

####

# Specify the location of the dorado bin and lib folders
DORADODIR=/home/stewarz2/various_programs/dorado/dorado-1.3.0-linux-x64

# Specify directory containing pod5 files for basecalling
POD5DIR=/work/ePGL/sequencing/dna/nanopore/citrus/NGS_TAL_Pete_041225/NGS_TAL_Pete_18Q053_041225/pod5

# Specify model/model-complex to use
MODEL=sup

# Specify output prefix and directory to write basecalled results to
OUTDIR=/work/ePGL/sequencing/dna/nanopore/citrus/NGS_TAL_Pete_041225/NGS_TAL_Pete_18Q053_041225/bam_sup
OUTPREFIX=18Q053

####

# STEP 0: Make sure LIB path is set and outdir exists
export LD_LIBRARY_PATH=${DORADODIR}/lib:$LD_LIBRARY_PATH
mkdir -p ${OUTDIR}

# STEP 0.5: Format model value for inclusion in output file name
MODELSUFFIX=$(echo ${MODEL} | tr , - | tr @ _)

# STEP 1: Run dorado
${DORADODIR}/bin/dorado basecaller ${MODEL} ${POD5DIR}/ > ${OUTDIR}/${OUTPREFIX}.${MODELSUFFIX}.bam
