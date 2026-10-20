#!/bin/bash -l
#PBS -N mgiraffe
#PBS -l walltime=24:00:00
#PBS -l mem=100G
#PBS -l ncpus=8
#PBS -J 1-X
#PBS -P <RPID>

cd $PBS_O_WORKDIR
set -euo pipefail

conda activate cactus
module load SAMtools/1.18 # get a working samtools, not conda's broken one

####

# Specify the giraffe reference location
REFDIR=/scratch/stewarz2/pangenomics/13Q019x18Q049/giraffe_reference
REFPREFIX=13Q019_18Q049_giraffe

# Specify reads dir
READSDIR=/scratch/stewarz2/pangenomics/13Q019x18Q049/reads
R1SUFFIX=.trimmed_1P.fq.gz
R2SUFFIX=.trimmed_2P.fq.gz

# Specify computational parameters
CPUS=8

####

TMPD=${PBS_O_WORKDIR}/tmp
mkdir -p ${TMPD}

# STEP 1: Find read file prefixes
declare -a READFILES
i=0
for f in ${READSDIR}/*${R1SUFFIX}; do
    READFILES[${i}]=$(echo "${f%%${R1SUFFIX}}");
    i=$((i+1));
done

# STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
FILEPREFIX=${READFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${FILEPREFIX})

# STEP 3: Map, surject to the reference paths, and add sample/read group
#   -N/-R set SM and RG so bcftools can tell samples apart
#   collate -> fixmate -> sort -> markdup is the standard duplicate-marking chain
vg giraffe -t ${CPUS} \
  -Z ${REFDIR}/${REFPREFIX}.d2.gbz \
  -m ${REFDIR}/${REFPREFIX}.d2.min \
  -d ${REFDIR}/${REFPREFIX}.d2.dist \
  -f ${FILEPREFIX}${R1SUFFIX} -f ${FILEPREFIX}${R2SUFFIX} \
  -N ${BASEPREFIX} -R ${BASEPREFIX} \
  -o BAM \
| samtools collate -@ ${CPUS} -O -u - ${TMPD}/${BASEPREFIX}.collate \
| samtools fixmate -@ ${CPUS} -m -u - - \
| samtools sort -@ ${CPUS} -m 2G -T ${TMPD}/${BASEPREFIX}.sort -u - \
| samtools markdup -@ ${CPUS} -T ${TMPD}/${BASEPREFIX}.md - ${BASEPREFIX}.giraffe.md.bam;

samtools index -@ ${CPUS} ${BASEPREFIX}.giraffe.md.bam

# STEP 4: Confirm BAM contig names match the reference FASTA exactly
diff <(samtools view -H ${BASEPREFIX}.giraffe.md.bam | awk '$1=="@SQ"{sub("SN:","",$2); print $2}' | sort) \
     <(cut -f1 ${REFDIR}/${REFPREFIX}.graphref.fa.fai | sort) > ${BASEPREFIX}.contig_check.txt \
  && echo "Contig names match reference" \
  || echo "WARNING: contig names differ from reference; see ${BASEPREFIX}.contig_check.txt";

# STEP 5: Quick mapping QC
samtools flagstat -@ ${CPUS} ${BASEPREFIX}.giraffe.md.bam > ${BASEPREFIX}.flagstat.txt
