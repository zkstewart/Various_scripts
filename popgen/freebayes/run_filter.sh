#!/bin/bash -l
#PBS -N vcffilter
#PBS -l walltime=00:30:00
#PBS -l mem=10G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

# Specify things needed to run this
GENOMEDIR=/home/stewarz2/flies/chapa_2021/genome
GENOME=btrys06_freeze2.rename.fasta

MAPDIR=/home/stewarz2/flies/chapa_2022/map


# > STEP 1: Get our file list
cd ${MAPDIR}
declare -a VCFFILES
i=0
for f in *.decomposed.vcf.gz; do
    VCFFILES[${i}]=$(echo "${f}");
    i=$((i+1));
done
cd $PBS_O_WORKDIR

# > STEP 2: Get our input files argument
SEPARATOR=" "
VCFTOOLS_ARG="$( printf "${SEPARATOR}%s" "${VCFFILES[@]}" )"

# > STEP 3: Merge individual VCFs
bcftools merge -Ov -o btrys06.merged.vcf ${VCFTOOLS_ARG}

# > STEP 3: Filter vcfs according to publication Pete emailed to pro_zac
MISSING=0.5 ## this means >50% of individuals need to have the site
MINQ=30 ## minimum SNP quality of 30, whatever that means
MAC=1 ## minor allele count must be >= 1
MINDEPTH=3 ## minimum depth of 3 for a genotype call
vcftools --vcf btrys06.merged.vcf --max-missing ${MISSING} --mac ${MAC} --minQ ${MINQ} --min-meanDP ${MINDEPTH} --remove-filtered-all --recode --recode-INFO-all --out btrys06.filtered.vcf
mv btrys06.filtered.vcf.recode.vcf btrys06.filtered.vcf

# > STEP 4: More custom filtering
POPMISSING=0.5 ## remove a site if each population does not have at least this much presence

