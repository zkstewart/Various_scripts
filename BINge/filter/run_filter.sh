#!/bin/bash -l
#PBS -N BINge_filter
#PBS -l walltime=08:00:00
#PBS -l mem=60G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

# Specify the location of the BINge directory
BINGEDIR=/home/stewarz2/scripts/BINge

# Specify the location of the BINge cluster file
CLUSTERDIR=/home/stewarz2/plant_group/ted/binge/BINge_ted_citrus
CLUSTERFILE=BINge_clustering_result.c37ec94651f78091ccc1b0d109c5af2bdf27ab7d115ae799a6c5a55623cc2c84.tsv

# Specify the location of the FASTA files and their format
FASTADIR=/home/stewarz2/plant_group/ted/binge/BINge_ted_citrus
SEQFORMAT=transcript # should be transcript, cds, or protein

# Specify the location of the BLAST file
BLASTDIR=/home/stewarz2/plant_group/ted/binge/blast
BLASTFILE=binge_mms2SEARCH_sorted.onlybest.m8

# Specify the location of the annotation file
ANNOTDIR=/home/stewarz2/plant_group/ted/binge
ANNOTFILE=citrus_annotations.gff3

# Specify the location of the salmon files
SALMONDIR=/home/stewarz2/plant_group/ted/binge/salmon
SALMONSUFFIX=.out

# Specify the output file name
OUTPUTFILE=BINge_ted_citrus.filtered.tsv

#####

# STEP 1: Get the salmon prefixes
declare -a SALMONNAMES
i=0
for dir in ${SALMONDIR}/*${SALMONSUFFIX}; do
        SALMONNAMES[${i}]=$(basename ${dir});
        i=$((i+1));
done

# STEP 2: Format the location of the salmon files
declare -a SALMONFILES
i=0
for name in ${SALMONNAMES[@]}; do
        SALMONFILES[${i}]=$(echo "${SALMONDIR}/${name}/quant.sf");
        i=$((i+1));
done

# STEP 3: Format salmon list for argument input
SEPARATOR=" "
SALMONFILE_ARG="$( printf "${SEPARATOR}%s" "${SALMONFILES[@]}" )"

# STEP 4: Run BINge_filter
python ${BINGEDIR}/BINge_filter.py -i ${CLUSTERDIR}/${CLUSTERFILE} \
	-f ${FASTADIR} -o ${OUTPUTFILE} \
	--annot ${ANNOTDIR}/${ANNOTFILE} \
	--blast ${BLASTDIR}/${BLASTFILE} \
	--salmon ${SALMONFILE_ARG} \
	-s ${SEQFORMAT} --require1x

