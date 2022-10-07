#!/bin/bash -l
#PBS -N merge_clair3
#PBS -l walltime=00:10:00
#PBS -l mem=15G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

#################################

# Specify the directory containing clair3 output VCFs
RESULTSDIR=main_results

# Specify the suffix that identifies fixed VCF files
SUFFIX=_fix.vcf.gz

# Specify the prefix for the output file
PREFIX=reticulata_clair3

#################################

# > STEP 1: Get our file list
declare -a VCFFILES
i=0
for f in ${RESULTSDIR}/*${SUFFIX}; do
    VCFFILES[${i}]=$(echo "${f}");
    i=$((i+1));
done

# > STEP 2: Get our input files argument
SEPARATOR=" "
VCFTOOLS_ARG="$( printf "${SEPARATOR}%s" "${VCFFILES[@]}" )"

# > STEP 3: Merge individual VCFs
bcftools merge -Ov -o ${PREFIX}.merged.vcf ${VCFTOOLS_ARG}

# > STEP 4: bgzip and index the VCF
bgzip -i ${PREFIX}.merged.vcf
bcftools index ${PREFIX}.merged.vcf.gz
gunzip -c ${PREFIX}.merged.vcf.gz > ${PREFIX}.merged.vcf
