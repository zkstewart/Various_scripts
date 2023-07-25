#!/bin/bash -l
#PBS -N BINgeQcav
#PBS -l walltime=60:00:00
#PBS -l mem=80G
#PBS -l ncpus=18
#PBS -W depend=afterok:4835131.pbs

cd $PBS_O_WORKDIR

# Specify the location of the BINge directory
BINGEDIR=/home/stewarz2/scripts/BINge

# Specify the location of the transcriptome file
FASTAFILES="/home/stewarz2/banana_group/qcav_dge/transcriptome/results/qcav_dge_okay-okalt.fasta /home/stewarz2/banana_group/annotations/Musa_acuminata_pahang_v4_main.trans /home/stewarz2/banana_group/annotations/Musa_balbisiana_main.trans"

# Specify the locations of the annotation files
ANNOTFILES="/home/stewarz2/banana_group/annotations/Musa_acuminata_pahang_v4.gff3 /home/stewarz2/banana_group/annotations/Musa_balbisiana.gff3"

# Specify the locations of the GMAP alignment files
GMAPFILES="/home/stewarz2/banana_group/qcav_dge/binge/gmap/qcav_to_A.6paths.gmap.gff3 /home/stewarz2/banana_group/qcav_dge/binge/gmap/qcav_to_B.6paths.gmap.gff3"

# Specify computational and behavioural parameters
MMSEQSDIR=/home/stewarz2/various_programs/mmseqs/bin
THREADS=18
OUTPUTFILE=BINge_qcav_clusters.tsv

#####

# Run BINge
python ${BINGEDIR}/BINge.py -i ${FASTAFILES} -ga ${ANNOTFILES} -gm ${GMAPFILES} -o ${OUTPUTFILE} \
    --threads ${THREADS} --clusterer mmseqs-cascade --mmseqs ${MMSEQSDIR}

