#!/bin/bash -l
#PBS -N sra
#PBS -l ncpus=1
#PBS -l walltime=16:00:00
#PBS -l mem=15G
#PBS -J 1-5

cd $PBS_O_WORKDIR

SCRATCHDIR=/scratch/stewarz2

####

TOOLKITDIR=/home/stewarz2/various_programs/sratoolkit.3.0.0-ubuntu64/bin
LISTFILE=datasets.txt

####

# STEP 0: Grab the SRA accession from the list file
ACC=$(cat ${LISTFILE} | head -n ${PBS_ARRAY_INDEX} | tail -n 1)

# STEP 1: Setup scratch location
mkdir -p ${SCRATCHDIR}/sratoolkit_${ACC}
cd ${SCRATCHDIR}/sratoolkit_${ACC}

# STEP 2: Download the file
$TOOLKITDIR/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files ${ACC}

# STEP 3: Gzip back to origin
gzip -c ${ACC}_1.fastq > ${PBS_O_WORKDIR}/${ACC}_1.fastq.gz
gzip -c ${ACC}_2.fastq > ${PBS_O_WORKDIR}/${ACC}_2.fastq.gz
