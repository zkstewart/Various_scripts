#!/bin/bash -l
#PBS -N qBINge_repr
#PBS -l walltime=06:00:00
#PBS -l mem=35G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

# Specify the location of the BINge directory
BINGEDIR=/home/stewarz2/scripts/BINge

# Specify the location of the BINge cluster file
CLUSTERDIR=/home/stewarz2/banana_group/qcav_dge/binge/filter
CLUSTERFILE=BINge_qcav_clusters.filtered.tsv # make sure this is the filtered file!

# Specify the location of the transcriptome file
FASTADIR1=/home/stewarz2/banana_group/qcav_dge/transcriptome/results
FASTAFILE1=qcav_dge_okay-okalt.aa

FASTADIR2=/home/stewarz2/banana_group/annotations
FASTAFILE2=Musa_acuminata_pahang_v4_main.aa

FASTADIR3=/home/stewarz2/banana_group/annotations
FASTAFILE3=Musa_balbisiana_main.aa

# Specify the location of the BLAST file
BLASTDIR=/home/stewarz2/banana_group/qcav_dge/binge/blast
BLASTFILE=binge_mms2SEARCH_sorted.onlybest.m8

# Specify the location of the annotation file
ANNOTDIR=/home/stewarz2/banana_group/metabolome/annotations
ANNOTFILE=Musa_combined.gff3

# Specify the location of the salmon files
SALMONDIR=/home/stewarz2/banana_group/qcav_dge/binge/salmon
SALMONNAMES="CL1.out CL2.out CL3.out CR1.out CR2.out CR3.out R_7L1.out R_7L2.out R_7L3.out R_7R1.out R_7R2.out R_7R3.out S_3L1.out S_3L2.out S_3L3.out S_3R1.out S_3R2.out S_3R3.out S_4L1.out S_4L2.out S_4L3.out S_4R1.out S_4R2.out S_4R3.out"

# Specify the output file name
OUTPUTFILE=qcav_representatives.BINge.filtered.aa


#####

# > STEP 1: Format the location of the salmon files
declare -a SALMONFILES
i=0
for name in $SALMONNAMES; do
        SALMONFILES[${i}]=$(echo "${SALMONDIR}/${name}/quant.sf");
        i=$((i+1));
done

# > STEP 2: Format salmon list for argument input
SEPARATOR=" "
SALMONFILE_ARG="$( printf "${SEPARATOR}%s" "${SALMONFILES[@]}" )"

# > STEP 3: Run BINge_representatives
python ${BINGEDIR}/BINge_representatives.py -i ${CLUSTERDIR}/${CLUSTERFILE} \
    -f ${FASTADIR1}/${FASTAFILE1} ${FASTADIR2}/${FASTAFILE2} ${FASTADIR3}/${FASTAFILE3} \
    -o ${OUTPUTFILE} --annot ${ANNOTDIR}/${ANNOTFILE} \
    --blast ${BLASTDIR}/${BLASTFILE} --salmon ${SALMONFILE_ARG}
