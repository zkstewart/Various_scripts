#!/bin/bash -l
#PBS -N star_index
#PBS -l walltime=02:00:00
#PBS -l mem=50G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

## MANUAL SETUP BELOW

# >> SETUP: Specify STAR executable location
STARDIR=/home/stewarz2/various_programs/STAR-2.7.10a/bin/Linux_x86_64_static

# >> SETUP: Specify genome location
GENDIR=/home/stewarz2/plant_group/leena/genome
GENFILE=NbLab360.genome.fasta

# >> SETUP: Specify gene annotation GTF location
## If you have a GFF3, convert it to GTF with agat_convert_sp_gff2gtf.pl
GTFDIR=/home/stewarz2/plant_group/leena/annotation
GTFFILE=NbLab360.v103.gtf

# >> SETUP: Specify computational resources
CPUS=1

## MANUAL SETUP END

## RUN PROGRAM
# STEP 1: Copy genome here and index it
cp ${GENDIR}/${GENFILE} .
${STARDIR}/STAR --runThreadN ${CPUS} \
	--runMode genomeGenerate \
	--genomeDir ${PBS_O_WORKDIR} \
	--genomeFastaFiles ${GENFILE} \
	--sjdbGTFfile ${GTFDIR}/${GTFFILE}
