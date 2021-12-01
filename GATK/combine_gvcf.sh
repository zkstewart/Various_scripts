#!/bin/bash -l
#PBS -N comb_02
#PBS -l walltime=04:00:00
#PBS -l mem=10G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

module load htslib/1.3-foss-2016a
module load bowtie2/2.2.9-foss-2016a
module load gatk/4.1.4.1-gcccore-8.3.0-java-11
module load atg/picard/2.2.2

export TILEDB_DISABLE_FILE_LOCKING=1

# MANUAL SETUP BELOW
## SETUP: Specify child directory to work within
CHILDDIR="02"

## SETUP: Specify gene annotation SAF
SAFDIR=/home/n8942188/plant_haplotypes/genes_of_interest
SAF=Alm_gois_final.saf

## SETUP: Specify genome file name and location
GENOMEDIR=/home/n8942188/plant_annotation/genomes/alm
GENOMEFILE=alm.fasta

## SETUP: Specify output directory
SPECIES_DIR=alm
# MANUAL SETUP END

## RUN PROGRAM
# STEP 1: Get our sample names
cd ${CHILDDIR}
declare -a SAMPLES
i=0
for d in *; do
    SAMPLES[${i}]=$(echo "${d}");
    i=$((i+1));
done

# STEP 2: Get our gene names
### Note: We assume all samples have the same gene files, violation of this assumption would be very bad
cd ${SAMPLES[0]}
declare -a GENES
i=0
for d in *; do 
    GENES[${i}]=$(echo "${d}");
    i=$((i+1));
done
cd ..

# STEP 3: Set up output directory structure
cd ${PBS_O_WORKDIR}
mkdir -p ${SPECIES_DIR}
cd ${SPECIES_DIR}

# STEP 4: Iterate through genes and run program
SEPARATOR=" -V "
for gene in "${GENES[@]}"; do
    # Find the contig this gene resides on
    CONTIG=$(echo "$(grep "${gene}" ${SAFDIR}/${SAF})" | awk '{print $2}')
    # Iterate through samples to get relevant VCF locations
    declare -a VCF_FILES=()
    i=0
    for sample in "${SAMPLES[@]}"; do
        FILE=${PBS_O_WORKDIR}/${CHILDDIR}/${sample}/${gene}/${gene}.vcf
        VCF_FILES[${i}]=$(echo "${FILE}");
        i=$((i+1));
    done
    # Format arguments to GATK
    V_ARGS="$( printf "${SEPARATOR}%s" "${VCF_FILES[@]}" )"
    # Call import on all VCFs
    gatk --java-options "-Xmx4g" GenomicsDBImport --genomicsdb-workspace-path ${gene}.db --intervals ${CONTIG} ${V_ARGS}
    # Run joint genotyping
    gatk --java-options "-Xmx4g" GenotypeGVCFs -R ${GENOMEDIR}/${GENOMEFILE} -V gendb://${gene}.db -new-qual -O ${gene}.genotypes.vcf
done
