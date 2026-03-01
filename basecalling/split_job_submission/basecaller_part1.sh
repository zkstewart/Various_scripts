#!/bin/bash -l
#PBS -N 662_p1
#PBS -l walltime=32:00:00
#PBS -l mem=15G
#PBS -l ncpus=1
#PBS -l ngpus=1
#PBS -l gpu_id=A100

cd $PBS_O_WORKDIR

module load zlib/1.3.1
export LD_LIBRARY_PATH=/mnt/weka/pkg/rhel94/AuthenticAMD-25/software/zlib/1.3.1/lib:${LD_LIBRARY_PATH}

####

# Specify the program locations
DORADODIR=/home/stewarz2/various_programs/dorado/dorado-1.3.0-linux-x64
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the prefix for output files
OUTPREFIX=NGS_662_

# Specify program behavioural values
MODEL=sup
MINREADLEN=5000

####

# STEP 0: Export and set all hard-coded variables
POD5DIR=pod5
BAMDIR=bam_sup
FQDIR=fastq
MODELSUFFIX=$(echo ${MODEL} | tr , - | tr @ _)
export LD_LIBRARY_PATH=${DORADODIR}/lib:${LD_LIBRARY_PATH}

# STEP 1: Run dorado
mkdir -p ${BAMDIR}
BAMRESULTFILE=${BAMDIR}/${OUTPREFIX}.${MODELSUFFIX}.bam
if [[ ! -f ${BAMRESULTFILE} ]]; then
    ${DORADODIR}/bin/dorado basecaller ${MODEL} ${POD5DIR}/ > ${BAMRESULTFILE};
fi
