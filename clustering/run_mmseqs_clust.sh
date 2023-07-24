#!/bin/bash -l
#PBS -N mms2Clust
#PBS -l walltime=78:00:00
#PBS -l mem=150G
#PBS -l ncpus=24

cd $PBS_O_WORKDIR

####

# SETUP: Specify the program locations
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts
MMS2DIR=/home/stewarz2/various_programs/mmseqs/bin

# SETUP: Specify the location of the input FASTA file
FASTADIR=/home/stewarz2/banana_group/metabolome/binge/unbinned
FASTA=tmp_BINge_unbinned_942b2d071be6a211ad63.fasta

# SETUP: Specify the output file name
OUTPUTFILE=qcav_unbinned_cascade_clusters.tsv

# SETUP: Specify clustering type
PROGRAM=mmseqs-cascade # mmseqs-cascade, mmseqs-linclust, or cd-hit

# SETUP: Specify program parameters
CPUS=24
IDENTITY=0.85 # mmseqs default == 0.9
COVERAGE=0.6 # mmseqs default == 0.8
EVALUE=1e-3 # mmseqs default == 1e-3
MODE=connected_component # set-cover, connected_component, or greedy
STEPS=3 # mmseqs default == 1; only relevant for cascade

####

# Run wrapper script
python ${VARSCRIPTDIR}/run_clustering.py -i ${FASTADIR}/${FASTA} -o ${OUTPUTFILE} -t ${CPUS} \
    -p ${PROGRAM} --mmseqs ${MMS2DIR} \
    --identity ${IDENTITY} --coverage ${COVERAGE} --evalue ${EVALUE} \
    --mode ${MODE} --steps ${STEPS}
