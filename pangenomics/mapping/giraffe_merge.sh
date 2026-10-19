#!/bin/bash -l
#PBS -N bcf_merge
#PBS -l walltime=06:00:00
#PBS -l mem=20G
#PBS -l ncpus=4

cd ${PBS_O_WORKDIR}

####

# Specify the prepared giraffe reference
REF=/scratch/stewarz2/pangenomics/13Q019x18Q049/giraffe_reference/13Q019_18Q049_giraffe.graphref.fa

# Specif the per-contig call files location
CONTIGSDIR=contig_calls

# Specify computational and behavioural params
CPUS=4

# Specify the output prefix
OUTPREFIX=13Q019_18Q049

####

# STEP 1: Concatenate per-contig VCFs in reference order
cut -f1 ${REF}.fai | tr '#' '_' | awk -v d=${INDIR} '{print d"/"$1".vcf.gz"}' > vcf_list.txt
bcftools concat --threads ${CPUS} -f vcf_list.txt -Oz -o ${OUTPREFIX}.raw.vcf.gz
bcftools index -t ${OUTPREFIX}.raw.vcf.gz

# STEP 2: Example site filters (tune using the distributions from bcftools stats)
bcftools filter --threads ${CPUS} \
  -s LowQual -e 'QUAL<30 || INFO/MQ<30' \
  ${OUTPREFIX}.raw.vcf.gz -Ou \
| bcftools filter --threads ${CPUS} \
  -S . -e 'FMT/DP<5 || FMT/GQ<20' \
  -Oz -o ${OUTPREFIX}.filtered.vcf.gz
bcftools index -t ${OUTPREFIX}.filtered.vcf.gz

# STEP 3: Summary stats
bcftools stats -s - ${OUTPREFIX}.filtered.vcf.gz > ${OUTPREFIX}.filtered.stats.txt
