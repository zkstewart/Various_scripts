#!/bin/bash -l
#PBS -N fgiraffe
#PBS -l walltime=06:00:00
#PBS -l mem=20G
#PBS -l ncpus=4
#PBS -P <RPID>

cd ${PBS_O_WORKDIR}
set -euo pipefail

####

# Specify filtering values
MISSING=0.5 ## this means >50% of individuals need to have the site
MINQ=30 ## minimum SNP quality of 30
MAC=1 ## minor allele count must be >= 1
MAF=0.05 ## minor allele frequency greater than or equal to 5%

# Specify the prefix of the VCF file
PREFIX=13Q019_18Q049

####

# STEP 1: Filter vcf
if [[ ! -f ${PREFIX}.wave.giraffe.filtered.vcf  ]]; then
  vcftools --gzvcf ${PREFIX}.wave.giraffe.vcf.gz --max-missing ${MISSING} --mac ${MAC} --minQ ${MINQ} --remove-filtered-all --recode --recode-INFO-all --maf ${MAF} --out ${PREFIX}.wave.giraffe.filtered.vcf;
  mv ${PREFIX}.wave.giraffe.filtered.vcf.recode.vcf ${PREFIX}.wave.giraffe.filtered.vcf;
  bgzip ${PREFIX}.wave.giraffe.filtered.vcf;
  tabix ${PREFIX}.wave.giraffe.filtered.vcf.gz;
fi
