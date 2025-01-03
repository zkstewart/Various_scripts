#!/bin/bash -l
#PBS -N bam2fq
#PBS -l walltime=03:30:00
#PBS -l mem=10G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify the location of the Various_scripts folder
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the location of the dorado bin and lib folders
DORADODIR=/home/stewarz2/various_programs/dorado/dorado-0.9.0-linux-x64

# Specify called BAM files from basecalling
BAMDIR=/work/ePGL/resequencing/nanopore/glauca_shearing_tests/Pete_13Q027_sample7_16092024/bam_sup
BAMPREFIX=glauca_shearing_13Q027
BAMSUFFIX=.sup.bam

# Specify output file details
FQDIR=fastq

# Specify behavioural parameters
MINREADLEN=5000

####

# STEP 1: Run samtools to convert BAM to FASTQ
mkdir -p ${FQDIR}
if [[ ! -f ${FQDIR}/${BAMPREFIX}.fastq ]]; then
    samtools fastq -T "*" ${BAMDIR}/${BAMPREFIX}${BAMSUFFIX} > ${FQDIR}/${BAMPREFIX}.fastq;
fi

# STEP 2: Filter FASTQ reads to minimum size cutoff
if [[ ! -f ${FQDIR}/${BAMPREFIX}.min${MINREADLEN}.fastq ]]; then
    python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f cullbelow -i ${FQDIR}/${BAMPREFIX}.fastq -n ${MINREADLEN} -o ${FQDIR}/${BAMPREFIX}.min${MINREADLEN}.fastq;
fi

# STEP 3: Figure out how many blocks we can chunk the FASTQ into for correction
NUMBLOCKS=$(${DORADODIR}/bin/dorado correct ${FQDIR}/${BAMPREFIX}.min${MINREADLEN}.fastq --compute-num-blocks)
echo "NUMBLOCKS = ${NUMBLOCKS}"
