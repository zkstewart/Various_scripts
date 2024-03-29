#!/bin/bash -l
#PBS -N dew_mango
#PBS -l walltime=240:00:00
#PBS -l mem=120G
#PBS -l ncpus=4

cd $PBS_O_WORKDIR

conda activate perl5
export PERL5LIB="/home/stewarz2/various_programs/DEW/PerlLib:$PERL5LIB"

# SETUP START
# >> SETUP: Secify Various_scripts dir
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# >> SETUP: Specify DEW location
DEWDIR=/home/stewarz2/various_programs/DEW

# >> SETUP: Specify DEW parameters
INSERT_SIZE=142

# >> SETUP: Specify prefix
SPECIES=mango_cultivars

# >> SETUP: Specify computational resources
CPUS=4

# >> SETUP: Gene models location
MODELDIR=/home/stewarz2/plant_group/mango_cultivars/annotation
MODELFILE=manindi_flc.nucl

# >> SETUP: Trimmomatic RNAseq reads dir
READSDIR=/home/stewarz2/plant_group/mango_cultivars/prepared_reads
SUFFIX=.fq.gz
BITTOTRIM=_

# SETUP END

#####

# RUN START
# >> STEP 1: Find RNAseq file prefixes
declare -a RNAFILES
i=0
for f in ${READSDIR}/*1${SUFFIX}; do
    RNAFILES[${i}]=$(echo "${f%%1${SUFFIX}}");
    i=$((i+1));
done

# >> STEP 2: Generate program arguments
SEPARATOR1="1${SUFFIX} "
READ1="$( printf "%s${SEPARATOR1}" "${RNAFILES[@]}" )"
READ1=${READ1::-1}

SEPARATOR2="2${SUFFIX} "
READ2="$( printf "%s${SEPARATOR2}" "${RNAFILES[@]}" )"
READ2=${READ2::-1}

SAMPLE_NAMES=""
for value in ${RNAFILES[@]}; do BASENAME=$(basename ${value} ${BITTOTRIM}); SAMPLE_NAMES+="${BASENAME},"; done
SAMPLE_NAMES=${SAMPLE_NAMES::-1}

# >> STEP 3: Pre-empt bugs present in DEW
## This is a bit finnicky, but we'll run with 2 samples, clean up the directory, fix the bug, then "resume" it
TMP_READ1="${RNAFILES[0]}1${SUFFIX} ${RNAFILES[1]}1${SUFFIX}"
TMP_READ2="${RNAFILES[0]}2${SUFFIX} ${RNAFILES[1]}2${SUFFIX}"
TMP_NAMES=""
for value in ${RNAFILES[@]:0:2}; do BASENAME=$(basename ${value} ${BITTOTRIM}); TMP_NAMES+="${BASENAME},"; done
TMP_NAMES=${TMP_NAMES::-1}
# >> 3.1: Dry run with two samples only
perl ${DEWDIR}/dew.pl -infile ${MODELDIR}/${MODELFILE} -format FASTQ -readset1 ${TMP_READ1} -readset2 ${TMP_READ2} -output ${SPECIES}_dew -threads ${CPUS} -kanga -remove_redund -no_pairwise -nographs -genomewide -readset_separation ${INSERT_SIZE} -only_alignments -sample_names ${TMP_NAMES}
# >> 3.2: Fix the bug in DEW
biokanga index -i ${PBS_O_WORKDIR}/${SPECIES}_dew_results/${SPECIES}_dew.toalign.noredundant.fsa_unaligned -o ${PBS_O_WORKDIR}/${SPECIES}_dew_results/${SPECIES}_dew.toalign.noredundant.fsa_unaligned.kangax -r ${PBS_O_WORKDIR}/${SPECIES}_dew_results/${SPECIES}_dew.toalign.noredundant.fsa_unaligned -t ${PBS_O_WORKDIR}/${SPECIES}_dew_results/${SPECIES}_dew.toalign.noredundant.fsa_unaligned
# >> 3.3: Clean up all other files
rm ${SPECIES}_dew_transcriptome.sqlite*
rm -r sort_tmp

cd ${SPECIES}_dew_results
rm ${SPECIES}_dew.toalign.noredundant.fsa_vs_*
rm -r edgeR
rm -r gene_coverage_plots
cd ..

# >> STEP 4: Run DEW alignment procedure
perl ${DEWDIR}/dew.pl -infile ${MODELDIR}/${MODELFILE} -format FASTQ -readset1 ${READ1} -readset2 ${READ2} -output ${SPECIES}_dew -threads ${CPUS} -kanga -remove_redund -no_pairwise -nographs -genomewide -readset_separation ${INSERT_SIZE} -only_alignments -sample_names ${SAMPLE_NAMES} -resume -over #-no_check

# >> STEP 5: Tabulate and summarise output read counts
python ${VARSCRIPTDIR}/DGE/mapping/merge_dew_stats.py -i ${SPECIES}_dew_results -o ${SPECIES}_dew_counts.tsv

