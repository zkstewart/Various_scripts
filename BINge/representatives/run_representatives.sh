#!/bin/bash -l
#PBS -N BINge_repr
#PBS -l walltime=06:00:00
#PBS -l mem=35G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

# Specify the location of the BINge directory
BINGEDIR=/home/stewarz2/scripts/BINge

# Specify the location of the BINge cluster file
CLUSTERDIR=/home/stewarz2/plant_group/ted/binge/filter
CLUSTERFILE=BINge_ted_citrus.filtered.tsv # make sure this is the filtered file!

# Specify the location of the FASTA files and their format
FASTADIR=/home/stewarz2/plant_group/ted/binge/BINge_ted_citrus

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
OUTPUTFILE=BINge_ted_citrus.BINge.filtered.trans

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

# STEP 4: Run BINge_representatives
python ${BINGEDIR}/BINge_representatives.py -i ${CLUSTERDIR}/${CLUSTERFILE} \
	-f ${FASTADIR} -o ${OUTPUTFILE} \
	--annot ${ANNOTDIR}/${ANNOTFILE} \
	--blast ${BLASTDIR}/${BLASTFILE} \
	--salmon ${SALMONFILE_ARG}
