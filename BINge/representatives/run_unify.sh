#!/bin/bash -l
#PBS -N BINge_unify
#PBS -l walltime=03:00:00
#PBS -l mem=25G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

# Specify the location of the BINge directory
BINGEDIR=/home/stewarz2/scripts/BINge

# Specify the location of the representatives FASTA
REPDIR=/home/stewarz2/banana_group/qcav_dge/binge/representatives
REPFILE=qcav_representatives.BINge.filtered.aa

# Specify the location of the file1 files
FASTA1DIR1=/home/stewarz2/banana_group/qcav_dge/transcriptome/results
FASTA1FILE1=qcav_dge_okay-okalt.aa

FASTA1DIR2=/home/stewarz2/banana_group/annotations
FASTA1FILE2=Musa_acuminata_pahang_v4_main.aa

FASTA1DIR3=/home/stewarz2/banana_group/annotations
FASTA1FILE3=Musa_balbisiana_main.aa

# Specify the location of the file2 files
FASTA2DIR1=/home/stewarz2/banana_group/qcav_dge/transcriptome/results
FASTA2FILE1=qcav_dge_okay-okalt.cds

FASTA2DIR2=/home/stewarz2/banana_group/annotations
FASTA2FILE2=Musa_acuminata_pahang_v4_main.nucl

FASTA2DIR3=/home/stewarz2/banana_group/annotations
FASTA2FILE3=Musa_balbisiana_main.nucl

# Specify the output file prefix
OUTPUTPREFIX=banana_representatives.BINge.filtered

#####

# Run unification
python ${BINGEDIR}/utilities/unify_representatives.py -r ${REPDIR}/${REPFILE} \
    -f1 ${FASTA1DIR1}/${FASTA1FILE1} ${FASTA1DIR2}/${FASTA1FILE2} ${FASTA1DIR3}/${FASTA1FILE3} \
    -f2 ${FASTA2DIR1}/${FASTA2FILE1} ${FASTA2DIR2}/${FASTA2FILE2} ${FASTA2DIR3}/${FASTA2FILE3} -o ${OUTPUTPREFIX}
