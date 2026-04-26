#!/bin/bash -l
#PBS -N gvcfpipe
#PBS -l walltime=12:00:00
#PBS -l mem=12G
#PBS -l ncpus=1
#PBS -J 1-X

cd $PBS_O_WORKDIR

####

# Specify the location of the genome FASTA
GENOMEDIR=/work/ePGL/genomes/citrus/clementina/jgi_v1/citrusgenomedb
GENOME=Cclementina_182_v1.fa

# Specify the location of the BAM files and their suffix
BAMDIR=/scratch/stewarz2/pollen_snps/NGS_532_James/mapping/jgi_v1
SUFFIX=.sorted.bam

####

# STEP 1: Get our file list
declare -a BAMFILES
i=0
for f in ${BAMDIR}/*${SUFFIX}; do
    BAMFILES[${i}]=$(echo "${f}");
    i=$((i+1));
done

# STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
BAMFILE=${BAMFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${BAMFILE} ${SUFFIX})

# STEP 3: Run variant calling pipeline
bcftools mpileup -Ou -f ${GENOMEDIR}/${GENOME} \
    -q 10 -Q 20 \
    --indels-cns --gvcf 0 \
    ${BAMFILE} | bcftools call -m --gvcf 0 -Oz -o ${BASEPREFIX}.gvcf.gz

# STEP 4: Normalise variant call positions
bcftools norm -f ${GENOMEDIR}/${GENOME} ${BASEPREFIX}.gvcf.gz > ${BASEPREFIX}.norm.gvcf.gz

# STEP 5: Index files
tabix ${BASEPREFIX}.norm.gvcf.gz
tabix -C ${BASEPREFIX}.norm.gvcf.gz
