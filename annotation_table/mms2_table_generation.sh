#!/bin/bash -l
#PBS -N refseqTable
#PBS -l walltime=06:00:00
#PBS -l mem=70G
#PBS -l ncpus=12

cd $PBS_O_WORKDIR

####

# Specify program locations
MMSEQDIR=/home/stewarz2/various_programs/mmseqs/bin
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify query FASTA location
QUERYDIR=/home/stewarz2/plant_group/avocado/flowering2021_remap/genome
QUERYFILE=GCA_029852735.1.protein.fasta

# Specify reference FASTA location
DB=RefSeq
DBDIR=/home/stewarz2/various_programs/refseq_db
DBFILE=refseq_plant.protein.faa

# Specify database annotation file locations
IM=/home/stewarz2/various_programs/uniref_db/idmapping_selected.tab
IO=/home/stewarz2/various_programs/uniref_db/go.obo

# Specify program parameters
CPUS=12
NUMHITS=15
NUM_ITERS=1

# Specify prefix for output files
OUTPREFIX=GCA_029852735.1_refseq

####

# STEP 0: Set PATH variable and tmp files directory
export PATH=${MMSEQDIR}/:$PATH
TMPDIR=${OUTPREFIX}_tmp

# STEP 1: Run MMseqs2
python ${VARSCRIPTDIR}/run_mmseqs2.py -q ${QUERYDIR}/${QUERYFILE} \
    -t ${DBDIR}/${DBFILE} \
    -o ${OUTPREFIX} \
    -m ${MMSEQDIR} -c ${CPUS} \
    -n ${NUM_ITERS} \
    --resume --blast_sort \
    --tmp ${TMPDIR} \
    -st auto

# STEP 2: Generate the annotation table
python ${VARSCRIPTDIR}/annotation_table/generate_annotation_table.py -rf ${QUERYDIR}/${QUERYFILE} \
    -bf ${DBDIR}/${DBFILE} -bo ${OUTPREFIX}_mms2SEARCH_sorted.m8 \
    -id ${IM} -io ${IO} \
    -o ${OUTPREFIX}.annotation.tsv \
    --numhits ${NUMHITS} \
    --db ${DB}
