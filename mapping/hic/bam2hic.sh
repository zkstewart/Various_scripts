#!/bin/bash -l
#PBS -N bam2hic
#PBS -l walltime=04:00:00
#PBS -l mem=20G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

conda activate hic # should contain pairtools

####

# Specify the juicer_tools .jar file
JAR=/home/stewarz2/various_programs/juicebox/juicer_tools.2.10.01.jar
THREEDDNA=/home/stewarz2/various_programs/3d-dna
JBSCRIPTS=/home/stewarz2/various_programs/juicebox/juicebox_scripts/juicebox_scripts

# Specify the reference genome contig sizes file
SIZES=/work/ePGL/genomes/citrus/superior/292/292.both_hap.sizes
GENOMEFASTA=/work/ePGL/genomes/citrus/superior/292/292.both_hap.fasta

# Specify the BAM file name
BAM=NGS743_292_S1.haphic.bam

# Specify output file prefix
PREFIX=292

####

# STEP 0: Create tmp dir
TMPDIR=${PBS_O_WORKDIR}/${PREFIX}_tmp
mkdir -p ${TMPDIR}

# STEP 1: Extract pairs from BAM, sort them, and save as a .pairs file
pairtools parse --chroms-path ${SIZES} ${BAM} | \
pairtools sort --tmpdir=${TMPDIR} -o ${PREFIX}.pairs

# STEP 2: Convert .pairs to .links format
grep -v '^#' ${PREFIX}.pairs |
awk '
BEGIN{OFS="\t"}
$2!="!" && $4!="!" {

    flag1=($6=="+")?0:16;
    flag2=($7=="+")?0:16;

    print flag1,$2,$3,0,
          flag2,$4,$5,
          1,1,"-","-",1,"-","-","-";
}' > ${PREFIX}.links

sort -k2,2 -k6,6 ${PREFIX}.links > ${PREFIX}.sorted.links

# STEP 3: Create .hic file from 3d-dna utilities
bash ${THREEDDNA}/visualize/run-assembly-visualizer.sh \
    -p false \
    ${PREFIX}.assembly \
    ${PREFIX}.sorted.links.txt

####

# Deprecated methdology below
## This approach filters out too much valuable information

# STEP 1: Create .assembly file
#python ${JBSCRIPTS}/makeAgpFromFasta.py ${GENOMEFASTA} ${PREFIX}.agp
#python ${JBSCRIPTS}/agp2assembly.py ${PREFIX}.agp ${PREFIX}.assembly

# STEP 2: Produce hic links
#matlock bam2 juicer ${BAM} ${PREFIX}.links.txt
#sort -k2,2 -k6,6 ${PREFIX}.links.txt > ${PREFIX}.sorted.links.txt

# STEP 3: Create .hic file
#bash ${THREEDDNA}/visualize/run-assembly-visualizer.sh \
#    -p false \
#    ${PREFIX}.assembly \
#    ${PREFIX}.sorted.links.txt

####

# Deprecated methdology below
## This approach does not render a .hic amenable to manual editing in juicebox

# STEP 0: Create tmp dir
#TMPDIR=${PBS_O_WORKDIR}/${PREFIX}_tmp
#mkdir -p ${TMPDIR}

# STEP 1: Extract pairs from BAM, sort them, and save as a .pairs file
#pairtools parse --chroms-path ${SIZES} ${BAM} | \
#pairtools sort --tmpdir=${TMPDIR} -o ${PREFIX}.pairs

# STEP 2: Run Juicer Tools for pairs -> hic conversion
#java -Xmx5g -jar ${JAR} pre ${PREFIX}.pairs ${PREFIX}.hic ${SIZES}
