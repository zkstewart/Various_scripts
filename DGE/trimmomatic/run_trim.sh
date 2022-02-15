#!/bin/bash -l
#PBS -N trim_alm
#PBS -l walltime=80:00:00
#PBS -l mem=50G
#PBS -l ncpus=1
#PBS -J 1-69

cd d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic

### MANUAL SETUP BELOW
## SETUP: Load modules
module load java/1.8.0_92

## SETUP: Specify trimmomatic location
TRIMDIR=D:\Bioinformatics\Protein_analysis\Trimmomatic-0.36
TRIMJAR=trimmomatic-0.36.jar

## SETUP: Specify file prefixes
SPECIES=alm

## SETUP: Specify RNAseq files
RNAFILES="d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/1_18-12-2019_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/2_18-12-2019_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/3_18-12-2019_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/5_18-12-2019_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/6_18-12-2019_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/7_18-12-2019_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/9_02-01-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/10_02-01-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/11_02-01-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/13_02-01-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/14_02-01-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/15_02-01-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/17_16-01-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/18_16-01-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/19_16-01-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/20_16-01-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/21_16-01-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/22_16-01-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/23_16-01-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/24_16-01-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/25_30-01-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/26_30-01-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/27_30-01-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/28_30-01-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/29_30-01-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/30_30-01-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/31_30-01-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/32_30-01-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/33_13-02-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/34_13-02-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/35_13-02-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/37_13-02-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/38_13-02-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/39_13-02-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/41_12-03-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/42_12-03-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/43_12-03-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/44_12-03-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/45_12-03-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/46_12-03-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/47_12-03-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/48_12-03-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/49_26-03-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/50_26-03-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/51_26-03-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/52_26-03-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/53_26-03-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/54_26-03-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/55_26-03-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/56_26-03-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/57_09-04-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/58_09-04-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/59_09-04-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/60_09-04-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/61_09-04-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/62_09-04-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/63_09-04-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/64_09-04-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/65_23-04-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/66_23-04-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/67_23-04-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/68_23-04-2020_Leaf d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/69_23-04-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/70_23-04-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/71_23-04-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/73_07-05-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/74_07-05-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/75_07-05-2020_Bud d:/Libraries/Documents/GitHub/Various_scripts/DGE/rnaseq_reads/76_07-05-2020_Bud"
ARRAY=($RNAFILES)
# Note: FILEPREFIXES assumes a file name like ${PREFIX}_1.fastq with paired ${PREFIX}_2.fastq
### MANUAL SETUP END

### AUTOMATIC SETUP START
## SETUP: Paired-end trimmomatic settings
COMMAND="ILLUMINACLIP:D:\Bioinformatics\Protein_analysis\Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"
### AUTOMATIC SETUP END

### RUN PROGRAM
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
FILE=${ARRAY[${ARRAY_INDEX}]}
BASENAME=$(basename ${FILE})

## STEP 1: Run Trimmomatic
java -jar $TRIMDIR/$TRIMJAR PE -threads 1 -trimlog $SPECIES.logfile ${FILE}_R1.fastq.gz ${FILE}_R2.fastq.gz -baseout d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/${BASENAME}.trimmed.fq.gz ${COMMAND}; done

## STEP 2: Unzip files
gunzip d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/${BASENAME}.trimmed_1P.fq.gz d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/${BASENAME}.trimmed_2P.fq; done
    