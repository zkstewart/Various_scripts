#!/bin/bash -l
#PBS -N pgiraffe
#PBS -l walltime=12:00:00
#PBS -l mem=60G
#PBS -l ncpus=8
#PBS -P <RPID>

cd $PBS_O_WORKDIR
set -euo pipefail

conda activate cactus
module load SAMtools/1.18 # get a working samtools, not conda's broken one

####

# Specify pangenome location
PANGENOMEDIR=/scratch/stewarz2/pangenomics/13Q019x18Q049/cactus_output
PANGENOMEPREFIX=13Q019_18Q049_pangenome

# Specify output location
REFDIR=/scratch/stewarz2/pangenomics/13Q019x18Q049/giraffe_reference
OUTPREFIX=13Q019_18Q049_giraffe

# Specify computational parameters
CPUS=8

####

REFSAMPLE=reference
GBZ=${REFDIR}/${OUTPREFIX}.d2.gbz
REF=${REFDIR}/${OUTPREFIX}.graphref.fa

# STEP 0: Set up the new reference location
mkdir -p ${REFDIR}
cp ${PANGENOMEDIR}/${PANGENOMEPREFIX}.d2.gbz ${GBZ}
echo "${PANGENOMEDIR}/${PANGENOMEPREFIX}.d2.gbz -> ${GBZ}" > ${REFDIR}/ORIGINAL_PANGENOME_LOCATION.txt

# STEP 1: List the reference-sample paths in the graph (sanity check)
vg paths -x ${GBZ} -L -S ${REFSAMPLE} > ${REFDIR}/ref_paths.txt

# STEP 2: Write those paths as FASTA and index it
vg paths -x ${GBZ} -F -S ${REFSAMPLE} > ${REF}
samtools faidx ${REF}

# STEP 3: Report the number of contigs (use this for the -J range in 03_bcftools_call.sh)
echo "Contigs: $(wc -l < ${REF}.fai)"

# STEP 4: Rebuild the vg index
vg version > ${REFDIR}/vg.version

vg index \
  -t ${CPUS} \
  -j ${REFDIR}/${OUTPREFIX}.d2.dist ${GBZ};

vg minimizer \
  -t ${CPUS} \
  -k 29 \
  -w 11 \
  -d ${REFDIR}/${OUTPREFIX}.d2.dist \
  -z ${REFDIR}/${OUTPREFIX}.d2.zipcodes \
  -o ${REFDIR}/${OUTPREFIX}.d2.min \
  ${GBZ};

# STEP 5: Reformat the pangenome for conformity to the later giraffe output
cut -f1 ${REF}.fai | awk '{n=$1; sub(/^reference#0#/,"",n); print n"\t"$1}' > ${REFDIR}/rename.txt

TEMP=${REFDIR}/giraffe_prep_tmp
mkdir -p ${TEMP}

bcftools annotate --rename-chrs ${REFDIR}/rename.txt \
  ${PANGENOMEDIR}/${PANGENOMEPREFIX}.wave.vcf.gz -Ou \
  | bcftools norm -f ${REF} -c e -m -any -Ou \
  | bcftools sort -T ${TEMP}/bcfsort -m 4G -Oz -o ${REFDIR}/${OUTPREFIX}.wave.renamed.norm.vcf.gz

bcftools index -t ${REFDIR}/${OUTPREFIX}.wave.renamed.norm.vcf.gz
