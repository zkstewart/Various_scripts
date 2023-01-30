#!/bin/bash -l
#PBS -N vcffil
#PBS -l walltime=04:00:00
#PBS -l mem=280G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

#################################

# Specify the location of the Various_scripts directory
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the location of the pop map file
POPSFILE=/home/stewarz2/flies/mitch/metadata/pops.txt

# Specify a prefix for output files
PREFIX=bactrocera

#################################


# > STEP 1: Filter vcfs according to publication Pete emailed to pro_zac
MISSING=0.5 ## this means >50% of individuals need to have the site
MINQ=30 ## minimum SNP quality of 30, whatever that means
MAC=1 ## minor allele count must be >= 1
MINDEPTH=3 ## minimum depth of 3 for a genotype call
MAF=0.05 ## minor allele frequency greater than or equal to 0.05
if [[ ! -f ${PREFIX}.filtered.vcf  ]]; then
    vcftools --vcf ${PREFIX}.decomposed.vcf --max-missing ${MISSING} --mac ${MAC} --minQ ${MINQ} --min-meanDP ${MINDEPTH} --remove-filtered-all --recode --recode-INFO-all --maf ${MAF} --out ${PREFIX}.filtered.vcf;
    mv ${PREFIX}.filtered.vcf.recode.vcf ${PREFIX}.filtered.vcf;
fi

# > STEP 2: Remove indels and keep only biallelic sites
if [[ ! -f ${PREFIX}.filtered.noindels.vcf  ]]; then
    bcftools view --max-alleles 2 --exclude-types indels -Ov -o ${PREFIX}.filtered.noindels.vcf ${PREFIX}.filtered.vcf
fi

# > STEP 3: More custom filtering
POPMISSING=0.5 ## remove a site if each population does not have at least this percentage of called genotypes
SAMPLEDP=5 ## set a site to be ambiguous if it does not have at least this much DP
if [[ ! -f ${PREFIX}.final.vcf  ]]; then
    python ${VARSCRIPTDIR}/popgen/VCF/filter_vcf.py -v ${PREFIX}.filtered.noindels.vcf -p ${POPSFILE} -o ${PREFIX}.final.vcf --mpp ${POPMISSING} --dp ${SAMPLEDP}
fi
