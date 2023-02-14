#!/bin/bash -l
#PBS -N vcfn_n
#PBS -l walltime=00:20:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -J 1-186

cd $PBS_O_WORKDIR

#################################

# Specify the location of the genome FASTA
GENOMEDIR=/home/stewarz2/flies/chapa_2021/genome
GENOME=btrys06_freeze2.rename.fasta

# Specify the mapping directory and file suffix
MAPDIR=/home/stewarz2/flies/chapa_2022/map
SUFFIX=.sorted.bam

#################################

# > STEP 1: Get our file list ## prefixes
cd ${MAPDIR}
declare -a BAMFILES
i=0
for f in *${SUFFIX}; do
    BAMFILES[${i}]=$(echo "${f}");
    i=$((i+1));
done
cd $PBS_O_WORKDIR

# > STEP 2: Get our array index
declare -i index
index=${PBS_ARRAY_INDEX}-1

# > STEP 3: Get our file for analysis
BAMFILE=${BAMFILES[${index}]}

# > STEP 4: Get our file prefix
PREFIX=${BAMFILE%%${SUFFIX}}

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
