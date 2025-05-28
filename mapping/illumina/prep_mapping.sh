#!/bin/bash -l

####

# Specify various scripts directory
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify location of BWA exe
BWA=/home/stewarz2/various_programs/bwa/bwa

# Specify genome file location
GENOMEDIR=/home/stewarz2/location/of/genome
GENOME=genome.fasta

# Specify fastq reads files location
FQDIR=/home/stewarz2/location/of/trimmed_reads

# Specify metadata csv location
CSV=/home/stewarz2/location/of/metadata/prep_metadata.csv

# Specify the column names
PREFIX=prefix # used to identify read file locations
ID=id # used to label your samples; VCF columns will contain this value as the header
SM=sample # this is a group label; if non-unique, reads in the same group will be variant called together; if you don't want that to happen, make this the same as ID

# Specify computational parameters
CPUS=4

####

# STEP 1: Generate shell script for job submission
python ${VARSCRIPTDIR}/mapping/illumina/illumina_map.py -b ${BWA} \
    -d ${FQDIR} -f ${GENOMEDIR}/${GENOME} -csv ${CSV} \
    --prefix ${PREFIX} --id ${ID} --sample ${SM} \
    --cpus ${CPUS}
