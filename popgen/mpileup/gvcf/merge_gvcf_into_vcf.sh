#!/bin/bash -l
#PBS -N mergegvcf
#PBS -l walltime=12:00:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify the location of the genome FASTA
GENOMEDIR=/work/ePGL/genomes/citrus/australasica/rainbow_v1
GENOME=rainbow_v1_hap2.fasta

# Specify the main VCF
VCF=/home/stewarz2/citrus/devindee/variants/CDS.filtered.vcf.gz

# Specify the location of any GVCF file(s) to merge
GVCF_DIR=/work/ePGL/sequencing/dna/illumina/citrus/NGS_532_James/gvcf/australasica/rainbow_v1

# Specify the output prefix
OUTPREFIX=flgenomics_commercial

####

# STEP 0: Locate normalised GVCF file(s)
find ${GVCF_DIR} -type f -name '*.norm.gvcf.gz' > gvcf_files.txt

# STEP 1: Merge GVCF(s) into main VCF file
bcftools merge --gvcf ${GENOMEDIR}/${GENOME} \
    -Oz -o ${OUTPREFIX}.gvcf.gz \
    --file-list gvcf_files.txt ${VCF}

# STEP 2: Filter resulting file to only list variant sites
bcftools view -c 1 ${OUTPREFIX}.gvcf.gz \
    -Oz -o ${OUTPREFIX}.gvcf_merged.vcf.gz
