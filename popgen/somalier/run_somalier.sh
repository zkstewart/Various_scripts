#!/bin/bash -l
#PBS -N somalier
#PBS -l walltime=04:00:00
#PBS -l mem=40G
#PBS -l ncpus=6

cd $PBS_O_WORKDIR

conda activate plink2

####

# Specify the VCF file location
VCF=/home/stewarz2/citrus/devindee/variants/merged.filtered.vcf.gz

# Specify the reference genome FASTA file location
GENOME=/home/stewarz2/citrus/devindee/genome/fl_hap1.fasta

# Specify the prefix for output files
PREFIX=flgenomics_commercial

####

# STEP 1: Add AF (allele frequency) tag into VCF
if [[ ! -f ${PREFIX}.step1.ok ]]; then
    bcftools +fill-tags ${VCF} \
        -Oz -o ${PREFIX}.af.vcf.gz --write-index \
        -- -t AN,AC,AF && touch ${PREFIX}.step1.ok;
fi;

# STEP 2: Run 'somalier find-sites'
if [[ ! -f ${PREFIX}.step2.ok ]]; then
    somalier find-sites \
        --output-vcf ${PREFIX}_sites.vcf.gz \
        --min-AN 2 ${PREFIX}.af.vcf.gz && touch ${PREFIX}.step2.ok;
fi;

# STEP 3: Run 'somalier extract'
if [[ ! -f ${PREFIX}.step3.ok ]]; then
    somalier extract -d ${PREFIX}_extracted/ \
        --sites ${PREFIX}_sites.vcf.gz \
        -f ${GENOME} \
        ${PREFIX}.af.vcf.gz && touch ${PREFIX}.step3.ok;
fi;

# STEP 4: Run 'somalier relate'
if [[ ! -f ${PREFIX}.step4.ok ]]; then
    somalier relate --infer \
        -o ${PREFIX}_somalier \
        ${PREFIX}_extracted/*.somalier && touch ${PREFIX}.step4.ok;
fi;
