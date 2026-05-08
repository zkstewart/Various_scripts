#!/bin/bash -l
#PBS -N star_index
#PBS -l walltime=02:00:00
#PBS -l mem=50G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

module load STAR/2.7.11b

####

# Specify genome file name
GENOME=alphonso_catas_2.1.fasta

# Specify gene annotation GTF location
## Convert GFF3 to GTF like: agat_convert_sp_gff2gtf.pl --gff input.gff3 -o output.gtf
GTF=alphonso_catas_2.1.gtf

# Specify computational resources
CPUS=1

# Specify output directory for the index
INDEXDIR=/work/ePGL/genomes/mango/indica/CATAS_Mindica_2.1/STAR_gtf_index

####

STAR --runThreadN ${CPUS} \
     --runMode genomeGenerate \
     --genomeDir ${INDEXDIR} \
     --genomeFastaFiles ${GENOME} \
     --sjdbGTFfile ${GTF}
