#!/bin/bash -l
#PBS -N pgcall
#PBS -l walltime=24:00:00
#PBS -l mem=20G
#PBS -l ncpus=4
#PBS -J 1-numContigs
#PBS -P <RPID>

## Get numContigs from: wc -l ${REF}.fai

cd ${PBS_O_WORKDIR}
set -euo pipefail

####

# Specify the prepared giraffe reference
REF=/scratch/stewarz2/pangenomics/13Q019x18Q049/giraffe_reference/13Q019_18Q049_giraffe.graphref.fa

# Specify the file listing all BAMs for joint calling
## Create this with "ls ${BAMDIR}/*${BAMSUFFIX} > ${BAMLIST}
BAMLIST=bamlist.txt

# Specify computational and behavioural params
CPUS=4
MINMAPQ=20
MINBASEQ=20
PLOIDY=2

# Specif the per-contig call files location
CONTIGSDIR=contig_calls

####

mkdir -p ${CONTIGSDIR}

# STEP 1: Get the contig for this subjob
REGION=$(sed -n "${PBS_ARRAY_INDEX}p" ${REF}.fai | cut -f1)
SAFENAME=$(echo "${REGION}" | tr '#' '_')

# STEP 2: Pileup -> multiallelic call -> left-align/normalise
bcftools mpileup --threads ${CPUS} \
  -f ${REF} -b ${BAMLIST} -r "${REGION}" \
  -q ${MINMAPQ} -Q ${MINBASEQ} \
  -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
  --max-depth 1000 -Ou \
| bcftools call --threads ${CPUS} -m -v --ploidy ${PLOIDY} -f GQ,GP -Ou \
| bcftools norm --threads ${CPUS} -f ${REF} -Oz -o ${CONTIGSDIR}/${SAFENAME}.vcf.gz

bcftools index -t ${CONTIGSDIR}/${SAFENAME}.vcf.gz
