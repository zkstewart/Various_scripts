#!/bin/bash -l

BAMDIRECTORY=/home/stewarz2/flies/chapa_2022/map
GENOMEFILE=/home/stewarz2/flies/chapa_2022/genome/btrys06_freeze2.rename.fasta
METADATAFILE=/home/stewarz2/flies/chapa_2022/metadata/popNums.txt
FREEBAYESEXE=/home/stewarz2/various_programs/freebayes/freebayes
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts
PREFIX=c2022
BAMSUFFIX=.sorted.bam
##
MODULES=modules.txt
SIZE=low

###########

#optional arguments:
#  -h, --help            show this help message and exit
#  -b BAMDIRECTORY       Input directory where BAM files reside
#  -g GENOMEFILE         Input genome fasta file
#  -m METADATAFILE       Input metadata text file
#  -f FREEBAYESEXE       Input the full location of the freebayes executable
#  -v VARSCRIPTDIR       Input the location of the Various_scripts directory
#  -p PREFIX             Specify the prefix for output files
#  -s BAMSUFFIX          Specify the suffix which identifies BAM files to be used
#  --modules MODULEFILE  Optionally provide a text file containing module load commands to embed within Freebayes scripts
#  --popMissing POPMISSING
#                        Optionally, during filtering, drop variant rows if each population does not have at least this proportion of individuals genotyped (default == 0.5)
#  --sampleDP SAMPLEDP   Optionally, during filtering, delete a sample's variant call if it's depth is lower than this value (default == 1)
#  --size {low,medium,high}
#                        Optionally, specify how much resources you think the jobs will need on a scale from low -> high (default == "medium")

###########

python ${VARSCRIPTDIR}/popgen/freebayes/individual/freebayes_persample_pipe.py -b ${BAMDIRECTORY} -g ${GENOMEFILE} -m ${METADATAFILE} -f ${FREEBAYESEXE} -v ${VARSCRIPTDIR} -p ${PREFIX} -s ${BAMSUFFIX} --modules ${MODULES} --size ${SIZE}

