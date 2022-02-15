#!/bin/bash -l
#PBS -N gvcf_02
#PBS -l walltime=04:00:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -J 1-17

cd $PBS_O_WORKDIR

module load htslib/1.3-foss-2016a
module load bowtie2/2.2.9-foss-2016a
module load gatk/4.1.4.1-gcccore-8.3.0-java-11
module load atg/picard/2.2.2

# MANUAL SETUP BELOW
## SETUP: Specify child directory to work within
CHILDDIR="02"

## SETUP: Specify genome file name and location
GENOMEDIR=/home/n8942188/plant_annotation/genomes/alm
GENOMEFILE=alm.fasta
# MANUAL SETUP END

# AUTOMATIC SETUP BELOW
## SETUP: Specify parent directory for RNAseq and mapping files
PARENT_OVERLAPDIR=/home/n8942188/plant_haplotypes/bt2_overlaps/${CHILDDIR}
PARENT_MAPDIR=/home/n8942188/plant_haplotypes/bt2_mapping/${CHILDDIR}

## SETUP: Specify the number of CPUS to use (always 1)
CPUS=1
# AUTOMATIC SETUP END

### NOTE: Make sure the genome file has a samtools index made with 'samtools faidx'

## RUN PROGRAM
# STEP 1: Get directory details to derive our samples
cd ${PARENT_OVERLAPDIR}
DIR=$(ls -lh | head -n $((PBS_ARRAY_INDEX+1)) | tail -n 1 | awk '{print $9}')

# STEP 2: Get our gene file names
cd ${PARENT_OVERLAPDIR}/${DIR}
declare -a GENES
i=0
for f in *_R1.fastq; do 
    GENES[${i}]=$(echo "${f%%_R1.fastq}");
    i=$((i+1));
done

# STEP 3: Set up directory structure
cd ${PBS_O_WORKDIR}
mkdir -p ${CHILDDIR}
cd ${CHILDDIR}
mkdir -p ${DIR}
cd ${DIR}

# STEP 4: Iterate through genes and run program
for gene in "${GENES[@]}"; do
    FILE1="${gene}_R1.fastq";
    FILE2="${gene}_R2.fastq";
    mkdir -p ${gene};
    cd ${gene};
    # Run Bowtie2 again for the selected reads
    bowtie2 -q -p ${CPUS} -N 1 -x ${PARENT_MAPDIR}/${GENOMEFILE}.bt2index -1 ${PARENT_OVERLAPDIR}/${DIR}/${FILE1} -2 ${PARENT_OVERLAPDIR}/${DIR}/${FILE2} -S ${gene}.sam
    # Convert SAM to BAM file
    samtools sort -m 5G -@ ${CPUS} -o ${gene}.sorted.bam -O bam ${gene}.sam;
    # Fix BAM read groups
    picard AddOrReplaceReadGroups I=${gene}.sorted.bam O=${gene}.sorted.fix.bam LB=lib1 PL=illumina PU=${GENOMEFILE} SM=${gene}
    # Index BAM file
    samtools index ${gene}.sorted.fix.bam
    # Call variants with GATK
    gatk --java-options "-Xmx4g" HaplotypeCaller -R ${GENOMEDIR}/${GENOMEFILE} -I ${gene}.sorted.fix.bam -O ${gene}.vcf --emit-ref-confidence GVCF;
    cd ..;
done
