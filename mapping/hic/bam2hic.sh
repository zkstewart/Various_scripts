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

# Specify the reference genome contig sizes file
SIZES=/work/ePGL/genomes/citrus/superior/292/292.both_hap.sizes

# Specify the BAM file name
BAM=NGS743_Minneola_S7.haphic.bam

# Specify output file prefix
PREFIX=292

####

# STEP 0: Create tmp dir
TMPDIR=${PBS_O_WORKDIR}/${PREFIX}_tmp
mkdir -p ${TMPDIR}

# STEP 1: Extract pairs from BAM, sort them, and save as a .pairs file
pairtools parse --chroms-path ${SIZES} ${BAM} | \
pairtools sort --tmpdir=${TMPDIR} -o ${PREFIX}.pairs

# STEP 2: Run Juicer Tools for pairs -> hic conversion
java -Xmx5g -jar ${JAR} pre ${PREFIX}.pairs ${PREFIX}.hic ${SIZES}
