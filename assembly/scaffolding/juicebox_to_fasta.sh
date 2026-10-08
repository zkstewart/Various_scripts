#!/bin/bash -l
#PBS -N jb2fasta
#PBS -l walltime=00:20:00
#PBS -l mem=15G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify the script repository locations
JBSCRIPTS=/home/stewarz2/various_programs/juicebox/juicebox_scripts/juicebox_scripts
ANNOTARIUM=/home/stewarz2/scripts/annotarium
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the input file locations
GENOMEFASTA=/work/ePGL/genomes/citrus/superior/292/292.both_hap.fasta
REVIEWASSEMBLY=292.review.assembly

# Specify the reference genome location
REFDIR=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/Pete_Bush_28_10_2024/reference_assembly
REFHAP1=Busk.hap1.fasta
REFHAP2=Busk.hap2.fasta

# Specify the prefix for output generation
PREFIX=292

####

# STEP 1: Convert review assembly back to fasta
python ${JBSCRIPTS}/juicebox_assembly_converter.py \
    -a ${REVIEWASSEMBLY} \
    -f ${GENOMEFASTA} \
    -p ${PREFIX} --simple_chr_names

# STEP 2: Filter out short fragments
python ${VARSCRIPTDIR}/fasta_handling_master_code.py \
    -f cullbelow \
    -i ${PREFIX}.fasta \
    -n 10000000 \
    -o ${PREFIX}.chrom.fasta

# STEP 3: Format a old:new IDs pair file
python ${ANNOTARIUM}/annotarium.py fasta ids \
    -i ${GENOMEFASTA} \
    -o ${PREFIX}.old.ids

python ${ANNOTARIUM}/annotarium.py fasta ids \
    -i ${PREFIX}.chrom.fasta \
    -o ${PREFIX}.new.ids

paste ${PREFIX}.new.ids ${PREFIX}.old.ids > ${PREFIX}.renaming.ids

# STEP 4: Rename the contigs
python ${ANNOTARIUM}/annotarium.py fasta rename \
    -i ${PREFIX}.chrom.fasta \
    -o ${PREFIX}.hic_edited.fasta \
    --sub ${PREFIX}.renaming.ids

# STEP 5: Separate out the two haplotypes again
python ${VARSCRIPTDIR}/fasta_handling_master_code.py \
    -f retrieveseqidwstring \
    -i ${PREFIX}.hic_edited.fasta \
    -s _hap1 \
    -o ${PREFIX}.hap1.fasta

python ${VARSCRIPTDIR}/fasta_handling_master_code.py \
    -f retrieveseqidwstring \
    -i ${PREFIX}.hic_edited.fasta \
    -s _hap2 \
    -o ${PREFIX}.hap2.fasta
