#!/bin/bash -l
#PBS -N vcf2plink
#PBS -l walltime=04:00:00
#PBS -l mem=40G
#PBS -l ncpus=6

cd $PBS_O_WORKDIR

conda activate plink2

####

# Specify the VCF file location
VCF=/home/stewarz2/citrus/devindee/variants/merged.filtered.vcf.gz

# Specify the FAM file location
FAM=/home/stewarz2/citrus/devindee/metadata/flgenomics_commercial.fam

# Specify the prefix for output files
PREFIX=flgenomics_commercial

# Specify behavioural parameters
MAXSNPMISS=0.1 # 10% missing data allowed
MAF=0.01 # 1% minor allele frequency minimum

####

# STEP 1: Run PLINK2 without setting --new-id-max-allele-len
PLINKERROR=${PREFIX}.tmpplink.error

if [[ ! -f ${PREFIX}.step1.ok ]]; then
    plink2 --vcf ${VCF} \
           --fam ${FAM} \
           --sort-vars --set-all-var-ids "@:#\$r,\$a" \
           --geno ${MAXSNPMISS} --maf ${MAF} \
           --make-pgen --out tmp_${PREFIX} 2> ${PLINKERROR} || touch ${PREFIX}.step1.ok;
fi;

# STEP 2: Identify the appropriate value to set for --new-id-max-allele-len
REGEX="has length ([0-9]+)\."
ERRORTEXT=$(cat ${PLINKERROR})

if [[ ${ERRORTEXT} =~ ${REGEX} ]]; then
    MAXALEN=$(echo "--new-id-max-allele-len ${BASH_REMATCH[1]} ");
else
    MAXALEN="";
fi;

# STEP 3: Generate PED representation of the data
if [[ ! -f ${PREFIX}.step3.ok ]]; then
    plink2 --vcf ${VCF} \
           --fam ${FAM} \
           --sort-vars --set-all-var-ids "@:#\$r,\$a" --rm-dup \
           ${MAXALEN} \
           --geno ${MAXSNPMISS} --maf ${MAF} \
           --make-pgen --out ${PREFIX} && touch ${PREFIX}.step3.ok;
fi;

# STEP 4: Generate BED representation of the data
if [[ ! -f ${PREFIX}.step4.ok ]]; then
    plink2 --vcf ${VCF} \
           --fam ${FAM} \
           --sort-vars --set-all-var-ids "@:#\$r,\$a" --rm-dup \
           ${MAXALEN} \
           --geno ${MAXSNPMISS} --maf ${MAF} \
           --make-bgen --out ${PREFIX} && touch ${PREFIX}.step4.ok;
fi;
