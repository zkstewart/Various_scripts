#!/bin/bash -l
#PBS -N vcffilter
#PBS -l walltime=08:00:00
#PBS -l mem=90G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

# Specify things needed to run this
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

GENOMEDIR=/home/stewarz2/flies/chapa_2021/genome
GENOME=btrys06_freeze2.rename.fasta

MAPDIR=/home/stewarz2/flies/chapa_2022/map

POPSFILE=/home/stewarz2/flies/chapa_2022/pops.txt

PREFIX=btrys06

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
if [[ ! -f ${PREFIX}.merged.vcf  ]]; then
    bcftools merge -Ov -o ${PREFIX}.merged.vcf ${VCFTOOLS_ARG}
fi

# > STEP 3: Filter vcfs according to publication Pete emailed to pro_zac
MISSING=0.5 ## this means >50% of individuals need to have the site
MINQ=30 ## minimum SNP quality of 30, whatever that means
MAC=1 ## minor allele count must be >= 1
MINDEPTH=3 ## minimum depth of 3 for a genotype call
MAF=0.05 ## minor allele frequency greater than or equal to 0.05
if [[ ! -f ${PREFIX}.filtered.vcf  ]]; then
    vcftools --vcf ${PREFIX}.merged.vcf --max-missing ${MISSING} --mac ${MAC} --minQ ${MINQ} --min-meanDP ${MINDEPTH} --remove-filtered-all --recode --recode-INFO-all --maf ${MAF} --out ${PREFIX}.filtered.vcf;
    mv ${PREFIX}.filtered.vcf.recode.vcf ${PREFIX}.filtered.vcf;
fi

# > STEP 4: Remove indels
if [[ ! -f ${PREFIX}.filtered.noindels.vcf  ]]; then
    bcftools view --exclude-types indels -Ov -o ${PREFIX}.filtered.noindels.vcf ${PREFIX}.filtered.vcf
fi

# > STEP 5: More custom filtering
POPMISSING=0.5 ## remove a site if each population does not have at least this much presence
if [[ ! -f ${PREFIX}.final.vcf  ]]; then
    python ${VARSCRIPTDIR}/popgen/VCF/filter_vcf.py -v ${PREFIX}.filtered.noindels.vcf -p ${POPSFILE} -o ${PREFIX}.final.vcf --mpp ${POPMISSING}
fi
if [[ ! -f ${PREFIX}.final.geno  ]]; then
    python ${VARSCRIPTDIR}/popgen/VCF/filter_vcf.py -v ${PREFIX}.filtered.noindels.vcf -p ${POPSFILE} --geno -o ${PREFIX}.final.geno --mpp ${POPMISSING}
fi
