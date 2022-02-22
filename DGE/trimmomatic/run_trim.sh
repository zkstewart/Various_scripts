#!/bin/bash -l
#PBS -N trim_kirra
#PBS -l walltime=24:00:00
#PBS -l mem=25G
#PBS -l ncpus=2
#PBS -J 1-4

cd $PBS_O_WORKDIR

### MANUAL SETUP BELOW
## SETUP: Load modules
module load java/1.8.0_92

## SETUP: Specify trimmomatic location
TRIMDIR=/home/stewarz2/various_programs/Trimmomatic-0.36
TRIMJAR=trimmomatic-0.36.jar

## SETUP: Specify file prefixes
SPECIES=btk

## SETUP: Specify RNAseq reads dir
READSDIR=/home/stewarz2/flies/kirra/reads/concat
SUFFIX=.fastq.gz

## SETUP: Specify computational parameters
CPUS=2
# SETUP END
### MANUAL SETUP END


COMMAND="ILLUMINACLIP:/home/stewarz2/various_programs/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"


#####

# RUN START
## STEP 1: Find RNAseq file prefixes
declare -a RNAFILES
i=0
for f in ${READSDIR}/*/*1${SUFFIX}; do
    RNAFILES[${i}]=$(echo "${f%%_R1${SUFFIX}}");
    i=$((i+1));
done

## STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
FILEPREFIX=${RNAFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${FILEPREFIX})

## STEP 1: Run Trimmomatic
java -jar $TRIMDIR/$TRIMJAR PE -threads ${CPUS} -trimlog ${SPECIES}.logfile ${FILEPREFIX}_R1${SUFFIX} ${FILE}_R2${SUFFIX} -baseout ${BASEPREFIX}.trimmed.fq.gz ${COMMAND}

## STEP 2: Unzip files
gunzip ${BASEPREFIX}.trimmed_1P.fq.gz ${BASEPREFIX}.trimmed_2P.fq
    
