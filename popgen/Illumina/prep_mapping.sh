#!/bin/bash -l

# START SCRIPT SETUP

# > START MANUAL SPECIFICATION
# >> Specify various scripts directory
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# >> Specify location of BWA exe
BWA=/home/stewarz2/various_programs/bwa/bwa

# >> Specify genome file location
GENOMEDIR=/home/stewarz2/flies/mitch/genome
GENOME=GCF_016617805.1_CSIRO_BtryS06_freeze2_genomic.fna

# >> Specify demultiplexed fastq file location
FQDIR=/home/stewarz2/flies/mitch/prepared_reads

# >> Specify metadata csv location
CSV=/home/stewarz2/flies/mitch/metadata/prep_metadata_concat_sp.csv

# >> Specify the column names
PREFIX=prefix
ID=id
SM=sm

# > END MANUAL SPECIFICATION

# END SCRIPT SETUP

# START SCRIPT RUN

# > STEP 1: Generate shell script for job submission
python ${VARSCRIPTDIR}/popgen/Illumina/illumina_map.py -d ${FQDIR} -f ${GENOMEDIR}/${GENOME} -csv ${CSV} -b ${BWA} --prefix ${PREFIX} --id ${ID} --sample ${SM} --cpus 4

# END SCRIPT RUN

