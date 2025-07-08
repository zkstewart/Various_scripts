#!/bin/bash -l
#PBS -N vcfFilter
#PBS -l walltime=12:00:00
#PBS -l mem=60G
#PBS -l ncpus=1

cd ${PBS_O_WORKDIR}

####

# Specify the prefix of the VCF .gz and for outputs
PREFIX=merged

####

# STEP 1: Filter vcf
#MISSING=0.5 ## this means >50% of individuals need to have the site
MINQ=30 ## minimum SNP quality of 30
#MAC=1 ## minor allele count must be >= 1
#MAF=0.05 ## minor allele frequency greater than or equal to 5%
if [[ ! -f ${PREFIX}.filtered.vcf.gz  ]]; then
    vcftools --gzvcf ${PREFIX}.vcf.gz --minQ ${MINQ} --remove-filtered-all --recode --recode-INFO-all --out ${PREFIX}.filtered.vcf;
    mv ${PREFIX}.filtered.vcf.recode.vcf ${PREFIX}.filtered.vcf;
    bgzip ${PREFIX}.filtered.vcf;
    tabix ${PREFIX}.filtered.vcf.gz;
fi

# STEP 2: (Optionally) remove indels and keep only biallelic sites
#if [[ ! -f ${PREFIX}.filtered.noindel.biallelic.vcf.gz  ]]; then
#    bcftools view --max-alleles 2 --exclude-types indels -Oz -o ${PREFIX}.filtered.noindel.biallelic.vcf.gz ${PREFIX}.filtered.vcf.gz;
#fi
