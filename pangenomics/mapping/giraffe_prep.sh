#!/bin/bash -l
#PBS -N pgiraffe
#PBS -l walltime=12:00:00
#PBS -l mem=20G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

conda activate cactus

####

# Specify pangenome location
PANGENOMEDIR=/scratch/stewarz2/pangenomics/13Q019x18Q049/cactus_output
PANGENOMEPREFIX=13Q019_18Q049_pangenome

# Specify output location
REFDIR=/scratch/stewarz2/pangenomics/13Q019x18Q049/giraffe_reference
OUTPREFIX=13Q019_18Q049_giraffe

####

REFSAMPLE=reference

mkdir -p ${REFDIR}
GBZ=${PANGENOMEDIR}/${PANGENOMEPREFIX}.d2.gbz
 
# STEP 1: List the reference-sample paths in the graph (sanity check)
vg paths -x ${GBZ} -L -S ${REFSAMPLE} > ${REFDIR}/ref_paths.txt
 
# STEP 2: Write those paths as FASTA and index it
vg paths -x ${GBZ} -F -S ${REFSAMPLE} > ${REFDIR}/${OUTPREFIX}.graphref.fa
samtools faidx ${REFDIR}/${OUTPREFIX}.graphref.fa
 
# STEP 3: Report the number of contigs (use this for the -J range in 03_bcftools_call.sh)
echo "Contigs: $(wc -l < ${REFDIR}/${OUTPREFIX}.graphref.fa.fai)"
