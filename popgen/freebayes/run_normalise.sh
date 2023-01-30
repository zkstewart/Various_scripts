#!/bin/bash -l
#PBS -N vcfn_n
#PBS -l walltime=24:00:00
#PBS -l mem=70G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

#################################

# Specify the location of the genome FASTA
GENOMEDIR=/home/stewarz2/flies/mitch/genome
GENOME=tryoni.fasta

# Specify the file prefix
PREFIX=bactrocera

#################################

# > STEP 1: Split multiallelic records to biallelic
if [[ ! -f ${PREFIX}.vcf.gz  ]]; then
    bgzip ${PREFIX}.vcf
fi
if [[ ! -f ${PREFIX}.vcf.gz.csi  ]]; then
    bcftools index ${PREFIX}.vcf.gz
fi
if [[ ! -f ${PREFIX}.split.vcf.gz  ]]; then
    bcftools norm -m- -Oz -o ${PREFIX}.split.vcf.gz -N ${PREFIX}.vcf.gz # with or without -N?
fi

# > STEP 2: Rejoin biallic sites into multiallelic sites
if [[ ! -f ${PREFIX}.rejoin.vcf.gz  ]]; then
    bcftools norm -m+ -Oz -o ${PREFIX}.rejoin.vcf.gz -N ${PREFIX}.split.vcf.gz # with or without -N
fi

# > STEP 3: Left-align and normalise everything
if [[ ! -f ${PREFIX}.normalised.vcf  ]]; then
    bcftools norm -f ${GENOMEDIR}/${GENOME} -Ov -o ${PREFIX}.normalised.vcf ${PREFIX}.rejoin.vcf.gz
fi

# > STEP 4: vt decompose SNPs
if [[ ! -f ${PREFIX}.decomposed.vcf  ]]; then
    vt decompose_blocksub ${PREFIX}.normalised.vcf > ${PREFIX}.decomposed.vcf
fi

# > STEP 5: Make file ready for next steps
if [[ ! -f ${PREFIX}.decomposed.vcf.gz  ]]; then
    bgzip -c ${PREFIX}.decomposed.vcf > ${PREFIX}.decomposed.vcf.gz
fi
if [[ ! -f ${PREFIX}.decomposed.vcf.gz.csi  ]]; then
    bcftools index ${PREFIX}.decomposed.vcf.gz
fi
