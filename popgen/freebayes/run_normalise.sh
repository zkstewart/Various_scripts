#!/bin/bash -l
#PBS -N vcfnorm
#PBS -l walltime=00:10:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -J 1-159

cd $PBS_O_WORKDIR

#################################

# Specify the location of the genome FASTA
GENOMEDIR=/home/stewarz2/daniel/genome
GENOME=btrys06_chr1.fasta

# Specify the location of the mapped BAM files
MAPDIR=/home/stewarz2/daniel/rnaseq/chr1_map

# Specify the suffix that identifies mapped BAM files
SUFFIX=.sorted.md.bam

#################################

# > STEP 1: Get our file list
declare -a BAMFILES
i=0
for f in ${MAPDIR}/*${SUFFIX}; do
    BAMFILES[${i}]=$(echo "${f}");
    i=$((i+1));
done

# > STEP 2: Get our array index
declare -i index
index=${PBS_ARRAY_INDEX}-1

# > STEP 3: Get our file for analysis
BAMFILE=${BAMFILES[${index}]}

# > STEP 4: Get our file prefix
PREFIX=$(basename ${BAMFILE} ${SUFFIX})

# > STEP 5: Split multiallelic records to biallelic
if [[ ! -f ${PREFIX}.vcf.gz  ]]; then
    bgzip ${PREFIX}.vcf
fi
if [[ ! -f ${PREFIX}.vcf.gz.csi  ]]; then
    bcftools index ${PREFIX}.vcf.gz
fi
if [[ ! -f ${PREFIX}.split.vcf.gz  ]]; then
    bcftools norm -m- -Oz -o ${PREFIX}.split.vcf.gz -N ${PREFIX}.vcf.gz # with or without -N?
fi

# > STEP 6: Rejoin biallic sites into multiallelic sites
if [[ ! -f ${PREFIX}.rejoin.vcf.gz  ]]; then
    bcftools norm -m+ -Oz -o ${PREFIX}.rejoin.vcf.gz -N ${PREFIX}.split.vcf.gz # with or without -N
fi

# > STEP 7: Left-align and normalise everything
if [[ ! -f ${PREFIX}.normalised.vcf  ]]; then
    bcftools norm -f ${GENOMEDIR}/${GENOME} -Ov -o ${PREFIX}.normalised.vcf ${PREFIX}.rejoin.vcf.gz
fi

# > STEP 8: vt decompose SNPs
if [[ ! -f ${PREFIX}.decomposed.vcf  ]]; then
    vt decompose_blocksub ${PREFIX}.normalised.vcf > ${PREFIX}.decomposed.vcf
fi

mkdir -p uncompressed
cp ${PREFIX}.decomposed.vcf uncompressed/${PREFIX}.decomposed.vcf

if [[ ! -f ${PREFIX}.decomposed.vcf.gz  ]]; then
    bgzip ${PREFIX}.decomposed.vcf
fi
if [[ ! -f ${PREFIX}.decomposed.vcf.gz.csi  ]]; then
    bcftools index ${PREFIX}.decomposed.vcf.gz
fi

