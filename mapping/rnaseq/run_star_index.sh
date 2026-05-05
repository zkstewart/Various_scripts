#!/bin/bash -l
#PBS -N star_index
#PBS -l walltime=02:00:00
#PBS -l mem=50G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify STAR executable location
STARDIR=/home/stewarz2/various_programs/STAR-2.7.11b/bin

# Specify genome file name
GENOME=/home/stewarz2/plant_group/leena/genome/NbLab360.genome.fasta

# Specify gene annotation GTF location
## If you have a GFF3, convert it to GTF with agat_convert_sp_gff2gtf.pl
GTF=/home/stewarz2/plant_group/leena/annotation/NbLab360.v103.gtf

# Specify computational resources
CPUS=1

# Specify output directory for the index
INDEXDIR=STAR_gtf_index

####

${STARDIR}/STAR --runThreadN ${CPUS} \
                --runMode genomeGenerate \
                --genomeDir ${INDEXDIR} \
                --genomeFastaFiles ${GENOME} \
                --sjdbGTFfile ${GTFDIR}/${GTFFILE}
