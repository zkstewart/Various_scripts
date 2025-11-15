#!/bin/bash -l
#PBS -N h1star
#PBS -l walltime=24:00:00
#PBS -l mem=60G
#PBS -l ncpus=12
#PBS -J 1-26

cd $PBS_O_WORKDIR

SCRATCHDIR=/scratch/stewarz2

####

# Specify STAR executable directory
STARDIR=/home/stewarz2/various_programs/STAR-2.7.10a/bin/Linux_x86_64_static

# Specify the location of genome FASTA files
GENOMEDIR=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/hap1

# Specify reads for mapping
LEFTREAD=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/transcriptome/normalised_reads/left.norm.fq
RIGHTREAD=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/transcriptome/normalised_reads/right.norm.fq

# Specify computational resources
CPUS=12
MEM=60

####

# STEP 1: Locate all target genome files
declare -a GENOMEFILES
i=0
for f in ${GENOMEDIR}/*_hap1.fasta; do
    GENOMEFILES[${i}]=$(echo "${f%%.fasta}");
    i=$((i+1));
done

# STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
GENOMEPREFIX=${GENOMEFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${GENOMEPREFIX})

# STEP 3: Set up a scratch dir for processing
cd ${SCRATCHDIR}
mkdir -p ${BASEPREFIX}_star_mapping
cd ${BASEPREFIX}_star_mapping

# STEP 4: Run STAR read mapping
${STARDIR}/STAR --runThreadN ${CPUS} \
    --genomeDir ${GENOMEDIR}/${BASEPREFIX}_star \
    --readFilesIn ${LEFTREAD} ${RIGHTREAD} \
    --twopassMode Basic

# STEP 5: Sort results
SAMTOOLSTHREADMEM=$(echo "$(printf "%.0f" $(echo "(${MEM}*0.50)/${CPUS}"|bc -l))")

samtools sort -m ${SAMTOOLSTHREADMEM}G \
    -@ ${CPUS} \
    -o ${BASEPREFIX}_star.sorted.bam \
    -O bam \
    Aligned.out.sam
samtools index ${BASEPREFIX}_star.sorted.bam

# STEP 6: Move files back out of scratch
mv ${BASEPREFIX}_star.sorted.bam $PBS_O_WORKDIR
mv ${BASEPREFIX}_star.sorted.bam.bai $PBS_O_WORKDIR
