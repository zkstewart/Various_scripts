#!/bin/bash -l
#PBS -N corset
#PBS -l walltime=80:00:00
#PBS -l mem=50G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

## MANUAL SETUP BELOW

# >> SETUP: Specify corset executable location
CORSETDIR=/home/stewarz2/various_programs/corset-1.09-linux64

# >> SETUP: Specify computational resources
CPUS=8

# >> SETUP: Specify salmon outputs location and 
MAPDIR=/home/stewarz2/anemones/cassie/mapping

# >> SETUP: Specify sample details specifically for corset
MAPRESULTS="An1_before_0h_S26.trimmed.out An2_before_0h_S27.trimmed.out An3_before_0h_S28.trimmed.out \
An1_before_72h_S32.trimmed.out An2_before_72h_S33.trimmed.out An3_before_72h_S34.trimmed.out \
An1_with_0h_S38.trimmed.out An2_with_0h_S39.trimmed.out An3_with_0h_S40.trimmed.out \
An1_with_72h_S44.trimmed.out An2_with_72h_S45.trimmed.out An3_with_72h_S46.trimmed.out"
##
GROUPING="1,1,1,2,2,2,3,3,3,4,4,4"
##
NAMES="An1_before_0h_S26 An2_before_0h_S27 An3_before_0h_S28 \
An1_before_72h_S32 An2_before_72h_S33 An3_before_72h_S34 \
An1_with_0h_S38 An2_with_0h_S39 An3_with_0h_S40 \
An1_with_72h_S44 An2_with_72h_S45 An3_with_72h_S46"

## MANUAL SETUP END

# STEP 1: Format input arguments
NAMES=$(echo ${NAMES} | tr " " ,)

SUFFIX="/aux_info/eq_classes.txt"
MAPRESULTS=(${MAPRESULTS})
MAPRESULTS=( "${MAPRESULTS[@]/%/${SUFFIX}}" )

SEPARATOR=" "
INPUTFILES="$( printf "%s${SEPARATOR}" "${MAPRESULTS[@]}" )"
INPUTFILES=${INPUTFILES::-1}

# STEP 2: Run corset
${CORSETDIR}/corset -g ${GROUPING} \
	-n ${NAMES} \
	-i salmon_eq_classes \
	-p corset_results \
	${INPUTFILES}
