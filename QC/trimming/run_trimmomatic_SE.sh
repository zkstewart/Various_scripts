#!/bin/bash -l
#PBS -N trim_rnaseq
#PBS -l walltime=24:00:00
#PBS -l mem=25G
#PBS -l ncpus=2
#PBS -J 1-12

cd $PBS_O_WORKDIR

### MANUAL SETUP BELOW
## SETUP: Load modules
module load java/1.8.0_92

## SETUP: Specify trimmomatic location
TRIMDIR=/home/stewarz2/various_programs/Trimmomatic-0.36
TRIMJAR=trimmomatic-0.36.jar

## SETUP: Specify file prefixes
SPECIES=leena

## SETUP: Specify RNAseq reads dir
READSDIR=/home/stewarz2/plant_group/leena/rnaseq_reads
SUFFIX=.fq

## SETUP: Specify computational parameters
CPUS=2
# SETUP END
### MANUAL SETUP END


COMMAND="ILLUMINACLIP:/home/stewarz2/various_programs/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"
## It is usually a better idea to get the precise adapter you used for sequencing and specify that file instead!

#####

# RUN START
## STEP 1: Find RNAseq file prefixes
declare -a RNAFILES
i=0
for f in ${READSDIR}/*${SUFFIX}; do
    RNAFILES[${i}]=$(echo "${f}");
    i=$((i+1));
done

## STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
SE_FILE=${RNAFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${SE_FILE} ${SUFFIX})

## STEP 1: Run Trimmomatic
java -jar $TRIMDIR/$TRIMJAR SE -threads ${CPUS} -trimlog ${SPECIES}.logfile ${SE_FILE} ${BASEPREFIX}.trimmed.fq.gz ${COMMAND}

## STEP 2: Unzip files
gunzip ${BASEPREFIX}.trimmed.fq.gz

