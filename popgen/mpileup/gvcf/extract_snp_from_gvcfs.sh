#!/bin/bash -l
#PBS -N extractsnp
#PBS -l walltime=00:30:00
#PBS -l mem=20G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify the location of the genome FASTA
GENOMEDIR=/work/ePGL/genomes/citrus/clementina/jgi_v1/citrusgenomedb
GENOME=Cclementina_182_v1.fa

# Specify the location to extract as VCF
POSITION="scaffold_8:6142000-6143000"

# Specify the output prefix
OUTPREFIX=pollen_snps

####

# STEP 0: Locate all normalised GVCF files nested in this location
find ~+ -type f -name '*.norm.gvcf.gz' > gvcf_files.txt

# STEP 1: Merge GVCF limited to SNP position range indicated
bcftools merge --file-list gvcf_files.txt \
    -r ${POSITION} \
    --gvcf ${GENOMEDIR}/${GENOME} \
    -Oz -o ${OUTPREFIX}.gvcf.gz

# STEP 2: Filter data to only list variant sites
bcftools view -c 1 ${OUTPREFIX}.gvcf.gz \
    -Oz -o ${OUTPREFIX}.vcf.gz
