#!/bin/bash -l
#PBS -N star_alm
#PBS -l walltime=80:00:00
#PBS -l mem=90G
#PBS -l ncpus=1
#PBS -W depend=afterok:
#PBS -J 1-69

cd d:\Libraries\Documents\GitHub\Various_scripts\DGE\star_map

## MANUAL SETUP BELOW
# SETUP: Specify STAR location
STARDIR=D:\Bioinformatics\Protein_analysis\STAR-2.7.6a\bin\Linux_x86_64

# SETUP: Specify prefix
SPECIES=alm

# SETUP: Specify computational resources
CPUS=1

# SETUP: Genome location
GENDIR=d:\Libraries\Documents\GitHub\Various_scripts\DGE\genome
GENFILE=alm.fasta

## SETUP: Specify RNAseq files
RNAFILES="d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/1_18-12-2019_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/2_18-12-2019_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/3_18-12-2019_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/5_18-12-2019_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/6_18-12-2019_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/7_18-12-2019_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/9_02-01-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/10_02-01-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/11_02-01-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/13_02-01-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/14_02-01-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/15_02-01-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/17_16-01-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/18_16-01-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/19_16-01-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/20_16-01-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/21_16-01-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/22_16-01-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/23_16-01-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/24_16-01-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/25_30-01-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/26_30-01-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/27_30-01-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/28_30-01-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/29_30-01-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/30_30-01-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/31_30-01-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/32_30-01-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/33_13-02-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/34_13-02-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/35_13-02-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/37_13-02-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/38_13-02-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/39_13-02-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/41_12-03-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/42_12-03-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/43_12-03-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/44_12-03-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/45_12-03-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/46_12-03-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/47_12-03-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/48_12-03-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/49_26-03-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/50_26-03-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/51_26-03-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/52_26-03-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/53_26-03-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/54_26-03-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/55_26-03-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/56_26-03-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/57_09-04-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/58_09-04-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/59_09-04-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/60_09-04-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/61_09-04-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/62_09-04-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/63_09-04-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/64_09-04-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/65_23-04-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/66_23-04-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/67_23-04-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/68_23-04-2020_Leaf.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/69_23-04-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/70_23-04-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/71_23-04-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/73_07-05-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/74_07-05-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/75_07-05-2020_Bud.trimmed_ d:\Libraries\Documents\GitHub\Various_scripts\DGE\trimmomatic/76_07-05-2020_Bud.trimmed_"
ARRAY=($RNAFILES)
## MANUAL SETUP END

## RUN PROGRAM
# STEP 1: Copy genome here. Need to do this since STAR can only tolerate 1 index per directory...
cp $GENDIR/$GENFILE .

# STEP 2: Generate index
$STARDIR/source/STAR --runThreadN $CPUS --runMode genomeGenerate --genomeDir d:\Libraries\Documents\GitHub\Various_scripts\DGE\star_map --genomeFastaFiles $GENFILE

# STEP 3: Run 2-pass procedure for each sample
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
FILE=${ARRAY[${ARRAY_INDEX}]}
BASENAME=$(basename ${FILE} .trimmed_)
mkdir $BASENAME
cd $BASENAME
$STARDIR/source/STAR  --runThreadN $CPUS --genomeDir d:\Libraries\Documents\GitHub\Various_scripts\DGE\star_map --readFilesIn ${FILE}1P.fq ${FILE}2P.fq --twopassMode Basic
cd ..
