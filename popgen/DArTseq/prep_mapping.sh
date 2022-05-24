#!/bin/bash -l

#PBS -N prep_map
#PBS -l ncpus=1
#PBS -l walltime=00:10:00
#PBS -l mem=5G

cd $PBS_O_WORKDIR

# START SCRIPT SETUP


####
DATA=Parental_Selected
####


# > START MANUAL SPECIFICATION
# >> Specify various scripts directory
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# >> Specify location of BWA exe
BWA=/home/stewarz2/various_programs/bwa/bwa

# >> Specify genome file location
GENOMEDIR=/home/stewarz2/flies/chapa_2021/genome
GENOME=btrys06_freeze2.rename.fasta

# >> Specify demultiplexed fastq file location
FQDIR=/home/stewarz2/flies/genome_based_2022/original/reads/${DATA}

# >> Specify metadata csv location
CSV=/home/stewarz2/flies/genome_based_2022/original/metadata/${DATA}/prep_metadata_parental.csv

# >> Specify the column names
ID=targetid
G="genotype location"

# > END MANUAL SPECIFICATION

# END SCRIPT SETUP

# START SCRIPT RUN

# > STEP 1: Generate shell script for job submission
python ${VARSCRIPTDIR}/popgen/DArTseq/dartseq_map.py -d ${FQDIR} -f ${GENOMEDIR}/${GENOME} -csv ${CSV} -id ${ID} -g ${G} -b ${BWA}

# END SCRIPT RUN

