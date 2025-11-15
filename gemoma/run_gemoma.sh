#!/bin/bash -l
#PBS -N h1gemoma
#PBS -l walltime=48:00:00
#PBS -l mem=40G
#PBS -l ncpus=14
#PBS -J 1-26

cd $PBS_O_WORKDIR

conda activate gemoma

SCRATCHDIR=/scratch/stewarz2

####

# Specify the location of the reference files
REFGENOME=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/references/concatenated/Busk.ref.bothhap.fasta
REFGFF3=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/references/concatenated/Busk.ref.bothhap.gff3

# Specify the directory of target files
TARGETDIR=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/hap1
TARGETSUFFIX=_hap1.fasta

# Specify the directory of RNAseq evidence mapping files
BAMSDIR=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/mapping/rna/hap1
BAMSUFFIX=_hap1_star.sorted.bam

# Specify haplotype
HAP=hap1

# Specify parameters
CPUS=14
MEM=30

####

# STEP 1: Locate all target genome files
declare -a TARGETFILES
i=0
for f in ${TARGETDIR}/*${TARGETSUFFIX}; do
    TARGETFILES[${i}]=$(echo "${f%%${TARGETSUFFIX}}");
    i=$((i+1));
done

# STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
TARGETPREFIX=${TARGETFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${TARGETPREFIX})

# STEP 3: Run GeMoMa pipeline on target genome
mkdir -p ${SCRATCHDIR}/${BASEPREFIX}_${HAP}_gemoma
cd ${SCRATCHDIR}/${BASEPREFIX}_${HAP}_gemoma
GeMoMa -Xmx${MEM}g GeMoMaPipeline t=${TARGETDIR}/${BASEPREFIX}${TARGETSUFFIX} \
    a=${REFGFF3} g=${REFGENOME} \
    r=MAPPED ERE.m=${BAMSDIR}/${BASEPREFIX}${BAMSUFFIX} \
    threads=${CPUS} AnnotationFinalizer.p=${BASEPREFIX}_${HAP}

# STEP 4: Move results back to origin
mv final_annotation.gff ${BASEPREFIX}_${HAP}.gemoma.gff3
mv predicted_proteins.fasta ${BASEPREFIX}_${HAP}.gemoma.proteins.fasta
mv protocol_GeMoMaPipeline.txt ${BASEPREFIX}_${HAP}.protocol_GeMoMaPipeline.txt
mv reference_gene_table.tabular ${BASEPREFIX}_${HAP}.reference_gene_table.tabular
