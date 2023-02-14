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

# > STEP 1: Drop sites which have no coverage
if [[ ! -f ${PREFIX}.dp.vcf  ]]; then
    vcffilter -f "DP > 0" ${PREFIX}.vcf > ${PREFIX}.dp.vcf
fi
if [[ ! -f ${PREFIX}.dp.vcf.gz  ]]; then
    bgzip -c ${PREFIX}.dp.vcf > ${PREFIX}.dp.vcf.gz
fi
if [[ ! -f ${PREFIX}.dp.vcf.gz.csi  ]]; then
    bcftools index ${PREFIX}.dp.vcf.gz
fi

# > STEP 2: Split multiallelic records to biallelic
if [[ ! -f ${PREFIX}.split.vcf.gz  ]]; then
    bcftools norm -m- -Oz -o ${PREFIX}.split.vcf.gz -N ${PREFIX}.dp.vcf.gz # with or without -N?
fi

# > STEP 3: Rejoin biallic sites into multiallelic sites
if [[ ! -f ${PREFIX}.rejoin.vcf.gz  ]]; then
    bcftools norm -m+ -Oz -o ${PREFIX}.rejoin.vcf.gz -N ${PREFIX}.split.vcf.gz # with or without -N
fi

# > STEP 4: Left-align and normalise everything
if [[ ! -f ${PREFIX}.normalised.vcf  ]]; then
    bcftools norm -f ${GENOMEDIR}/${GENOME} -Ov -o ${PREFIX}.normalised.vcf ${PREFIX}.rejoin.vcf.gz
fi

# > STEP 5: vt decompose SNPs
if [[ ! -f ${PREFIX}.decomposed.vcf  ]]; then
    vt decompose_blocksub ${PREFIX}.normalised.vcf > ${PREFIX}.decomposed.vcf
fi

# > STEP 6: Make file ready for next steps
if [[ ! -f ${PREFIX}.decomposed.vcf.gz  ]]; then
    bgzip -c ${PREFIX}.decomposed.vcf > ${PREFIX}.decomposed.vcf.gz
fi
if [[ ! -f ${PREFIX}.decomposed.vcf.gz.csi  ]]; then
    bcftools index ${PREFIX}.decomposed.vcf.gz
fi
