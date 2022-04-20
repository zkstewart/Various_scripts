#!/bin/bash -l
#PBS -N vcfnorm
#PBS -l walltime=00:10:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -J 1-186

cd $PBS_O_WORKDIR

# Module loads for vt
module unload gcc/4.9.3-2.25
module unload gcccore/4.9.3
module unload binutils/2.25-gcccore-4.9.3
module unload zlib/1.2.8-foss-2016a
module load gcc/10.3.0

module unload libxml2/2.9.3-foss-2016a
module load libxml2/2.9.10-gcccore-10.3.0

# Specify things needed to run this
GENOMEDIR=/home/stewarz2/flies/chapa_2021/genome
GENOME=btrys06_freeze2.rename.fasta

MAPDIR=/home/stewarz2/flies/chapa_2022/map


# > STEP 1: Get our file list ## prefixes
cd ${MAPDIR}
declare -a BAMFILES
i=0
for f in *.sorted.bam; do
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
PREFIX=${BAMFILE%%.sorted.bam}

# > STEP 5: Split multiallelic records to biallelic
bgzip ${PREFIX}.vcf
bcftools index ${PREFIX}.vcf.gz
bcftools norm -m- -Oz -o ${PREFIX}.split.vcf.gz -N ${PREFIX}.vcf.gz # with or without -N?

# > STEP 6: Rejoin biallic sites into multiallelic sites
bcftools norm -m+ -Oz -o ${PREFIX}.rejoin.vcf.gz -N ${PREFIX}.split.vcf.gz # with or without -N

# > STEP 7: Left-align and normalise everything
bcftools norm -f ${GENOMEDIR}/${GENOME} -Ov -o ${PREFIX}.normalised.vcf ${PREFIX}.rejoin.vcf.gz

# > STEP 8: vt decompose SNPs
vt decompose_blocksub ${PREFIX}.normalised.vcf > ${PREFIX}.decomposed.vcf
cp ${PREFIX}.decomposed.vcf uncompressed/${PREFIX}.decomposed.vcf
bgzip ${PREFIX}.decomposed.vcf
bcftools index ${PREFIX}.decomposed.vcf.gz
