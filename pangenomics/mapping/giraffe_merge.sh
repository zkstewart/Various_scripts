#!/bin/bash -l
#PBS -N mgiraffe
#PBS -l walltime=06:00:00
#PBS -l mem=20G
#PBS -l ncpus=4
#PBS -P <RPID>

cd ${PBS_O_WORKDIR}
set -euo pipefail

####

# Specify the prepared giraffe reference location
REFDIR=/scratch/stewarz2/pangenomics/13Q019x18Q049/giraffe_reference
REFPREFIX=13Q019_18Q049_giraffe

# Specif the per-contig call files location
CONTIGSDIR=contig_calls

# Specify computational and behavioural params
CPUS=4

# Specify the output prefix
OUTPREFIX=13Q019_18Q049

####

TEMP=giraffe_merge_tmp
mkdir -p ${TEMP}

# STEP 1: Concatenate per-contig VCFs in reference order
cut -f1 ${REFDIR}/${REFPREFIX}.graphref.fa.fai | tr '#' '_' | awk -v d=${CONTIGSDIR} '{print d"/"$1".vcf.gz"}' > vcf_list.txt
bcftools concat --threads ${CPUS} -f vcf_list.txt -Oz -o ${OUTPREFIX}.raw.vcf.gz
bcftools index -t ${OUTPREFIX}.raw.vcf.gz

# STEP 2: Normalise the merged VCF analogously to how the giraffe reference was handled
bcftools norm -f ${REFDIR}/${REFPREFIX}.graphref.fa -c e -m -any -Ou ${OUTPREFIX}.raw.vcf.gz \
  | bcftools sort -T ${TEMP}/bcfsort -m 4G -Oz -o ${OUTPREFIX}.raw.normalised.vcf.gz
  bcftools index -t ${OUTPREFIX}.raw.normalised.vcf.gz

# STEP 3: Extract overlapping sites between the pangenome VCF and this giraffe VCF
bcftools isec -n=2 -w1 \
  -Oz -o ${TEMP}/${OUTPREFIX}.isec.vcf.gz \
  ${OUTPREFIX}.raw.normalised.vcf.gz ${REFDIR}/${REFPREFIX}.wave.renamed.norm.vcf.gz
bcftools index -t ${TEMP}/${OUTPREFIX}.isec.vcf.gz

# STEP 4: Merge the common sites
bcftools merge -m none \
  -Oz -o ${OUTPREFIX}.wave.giraffe.vcf.gz --write-index \
  ${TEMP}/${OUTPREFIX}.isec.vcf.gz ${REFDIR}/${REFPREFIX}.wave.renamed.norm.vcf.gz
bcftools index -t ${OUTPREFIX}.wave.giraffe.vcf.gz

# STEP 4: Summary stats
bcftools stats -s - ${OUTPREFIX}.raw.normalised.vcf.gz > ${OUTPREFIX}.raw.normalised.stats.txt
bcftools stats -s - ${OUTPREFIX}.wave.giraffe.vcf.gz > ${OUTPREFIX}.wave.giraffe.stats.txt

## Note for future development
# It is at this point that I realise we have a minor problem. We are not using GVCF
# so we cannot tell missing sites are invariant or non-existing (no coverage).
# In the future, I should run giraffe_call something like below:
# > bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' $W | bgzip > wave.alleles.tsv.gz
# > tabix -s1 -b2 -e2 wave.alleles.tsv.gz
# > bcftools mpileup -f $REF -b bams.txt -T wave.alleles.tsv.gz \
#  -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 -Ou \
# | bcftools call -m -C alleles -T wave.alleles.tsv.gz -f GQ -Oz -o giraffe.at_wave_alleles.vcf.gz
# Then, here I will not use bcftools isec but instead use bcftools merge
