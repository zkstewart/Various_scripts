#!/bin/bash -l
#PBS -N mms2_ppole
#PBS -l walltime=46:00:00
#PBS -l mem=750G
#PBS -l ncpus=32

cd $PBS_O_WORKDIR

####

# Specify program locations
MMSEQDIR=/home/stewarz2/various_programs/mmseqs/bin
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify query FASTA location
QUERYDIR=/home/stewarz2/citrus/andrew_miles/powerpole/genome
QUERYFILE=TEMP_hp.OGS2.protein.fasta

# Specify target FASTA location
DBDIR=/home/stewarz2/various_programs/uniref_db
DBFILE=uniref90.fasta

# Specify output details
OUTPREFIX=powerpole_hp

# Specify tmp directory
TMPDIR=hp_uniref

# Specify program parameters
CPUS=32
NUM_ITERS=1
SENS=5.7

####

# Step 1: Add MMSEQS to path
export PATH=${MMSEQDIR}/:$PATH

# Step 2: Run MMseqs2
python ${VARSCRIPTDIR}/run_mmseqs2.py -q ${QUERYDIR}/${QUERYFILE} \
    -t ${DBDIR}/${DBFILE} \
    -o ${OUTPREFIX} \
    -m ${MMSEQDIR} \
    -c ${CPUS} \
    -n ${NUM_ITERS} -s ${SENS} \
    --resume \
    --blast_sort \
    --tmpDir ${TMPDIR}
