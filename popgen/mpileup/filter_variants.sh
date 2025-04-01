#!/bin/bash -l
#PBS -N filtervcf
#PBS -l walltime=12:00:00
#PBS -l mem=40G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify filtering values
MISSING=0.5 ## this means >50% of individuals need to have the site
MINQ=30 ## minimum SNP quality of 30
MAC=1 ## minor allele count must be >= 1
MAF=0.05 ## minor allele frequency greater than or equal to 5%

# Specify the prefix of the VCF file
PREFIX=merged

####

# STEP 1: Filter vcf
if [[ ! -f ${PREFIX}.filtered.vcf  ]]; then
    vcftools --gzvcf ${PREFIX}.vcf.gz --max-missing ${MISSING} --mac ${MAC} --minQ ${MINQ} --remove-filtered-all --recode --recode-INFO-all --maf ${MAF} --out ${PREFIX}.filtered.vcf;
    mv ${PREFIX}.filtered.vcf.recode.vcf ${PREFIX}.filtered.vcf;
fi

# STEP 2: Keep only biallelic sites
if [[ ! -f ${PREFIX}.filtered.biallelic.vcf  ]]; then
    bcftools view --max-alleles 2 -Ov -o ${PREFIX}.filtered.biallelic.vcf ${PREFIX}.filtered.vcf
fi

# STEP 3: Remove indels
if [[ ! -f ${PREFIX}.filtered.biallelic.noindels.vcf  ]]; then
    bcftools view --exclude-types indels -Ov -o ${PREFIX}.filtered.biallelic.noindels.vcf ${PREFIX}.filtered.biallelic.vcf
fi
