#!/bin/bash -l
#PBS -N BINgecount
#PBS -l walltime=01:00:00
#PBS -l mem=20G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

# Specify the location of the BINge directory
BINGEDIR=/home/stewarz2/scripts/BINge

# Specify the location of the BINge cluster file
CLUSTERDIR=/home/stewarz2/banana_group/qcav_dge/binge/filter
CLUSTERFILE=BINge_qcav_clusters.filtered.tsv

# Specify the location of the salmon files
SALMONDIR=/home/stewarz2/banana_group/qcav_dge/binge/salmon
SALMONNAMES="CL1.out CL2.out CL3.out CR1.out CR2.out CR3.out R_7L1.out R_7L2.out R_7L3.out R_7R1.out R_7R2.out R_7R3.out S_3L1.out S_3L2.out S_3L3.out S_3R1.out S_3R2.out S_3R3.out S_4L1.out S_4L2.out S_4L3.out S_4R1.out S_4R2.out S_4R3.out"


# Specify the output file name
OUTPUTFILE=qcav_counts.BINge.filtered.onlyBinned


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

# > STEP 2: Run BINge_counter
python ${BINGEDIR}/BINge_counter.py -i ${CLUSTERDIR}/${CLUSTERFILE} \
    -s ${SALMONFILE_ARG} -n ${SALMONNAMES} \
    -o ${OUTPUTFILE} --only_binned
