#!/bin/bash -l
#PBS -N demux
#PBS -l walltime=06:00:00
#PBS -l mem=35G
#PBS -l ncpus=4
#PBS -J 1-X

cd $PBS_O_WORKDIR

conda activate twist # should have bioconda::fgbio installed

####

# Specify Various_scripts location
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify reads dir
READSDIR=/work/ePGL/sequencing/dna/illumina/mango/COY19345_Will/prepared_reads
R1SUFFIX=_1.fq.gz
R2SUFFIX=_2.fq.gz

####

# STEP 0: Locate the well_barcodes.csv file
BARCODES=${VARSCRIPTDIR}/demultiplexing/twist/well_barcodes.csv

# STEP 1: Find file prefixes
declare -a RNAFILES
i=0
for f in ${READSDIR}/*${R1SUFFIX}; do
    RNAFILES[${i}]=$(echo "${f%%${R1SUFFIX}}");
    i=$((i+1));
done

# STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
FILEPREFIX=${RNAFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${FILEPREFIX})

# STEP 3: Run demultiplexing
if [[ ! -f ${BASEPREFIX}.ok ]]; then
    fgbio DemuxFastqs \
        --metadata ${BARCODES} \
        --min-mismatch-delta 2 \
        --max-mismatches 1 \
        --read-structures 6B2S+T 6B2S+T \
        --inputs ${FILEPREFIX}${R1SUFFIX} ${FILEPREFIX}${R2SUFFIX} \
        --output ${BASEPREFIX} \
        --output-type Fastq \
        --metrics ${BASEPREFIX}.tsv && touch ${BASEPREFIX}.ok;
fi;
