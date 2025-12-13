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

# Specify the program locations
DORADODIR=/home/stewarz2/various_programs/dorado/dorado-1.3.0-linux-x64
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the prefix for output files
OUTPREFIX=18Q053

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

# STEP 2: Run samtools to convert BAM to FASTQ
mkdir -p ${FQDIR}
FQRESULTFILE=${FQDIR}/${OUTPREFIX}.fastq
if [[ ! -f ${FQRESULTFILE} ]]; then
    samtools fastq -T "*" ${BAMRESULTFILE} > ${FQRESULTFILE};
fi

# STEP 3: Filter FASTQ reads to minimum size cutoff
FQFILTERFILE=${FQDIR}/${OUTPREFIX}.min${MINREADLEN}.fastq
if [[ ! -f ${FQDIR}/${OUTPREFIX}.min${MINREADLEN}.fastq ]]; then
    python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f cullbelow -i ${FQRESULTFILE} -n ${MINREADLEN} -o ${FQFILTERFILE};
fi

# STEP 4: Figure out how many blocks we can chunk the FASTQ into for correction
NUMBLOCKS=$(${DORADODIR}/bin/dorado correct ${FQFILTERFILE} --compute-num-blocks)
echo "NUMBLOCKS = ${NUMBLOCKS}"
echo "${NUMBLOCKS}" > ${FQDIR}/numblocks.txt
