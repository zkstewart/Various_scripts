#!/bin/bash -l
#PBS -N mms2_
#PBS -l walltime=120:00:00
#PBS -l mem=750G
#PBS -l ncpus=8

cd $PBS_O_WORKDIR


#################################

# Specify program locations
MMSEQDIR=/home/stewarz2/various_programs/MMseqs2/mmseqs/bin
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify query FASTA location
QUERYDIR=/home/stewarz2/citrus/annotation
QUERYFILE=GJ.v1.0.gene.model.aa

# Specify target FASTA location
DBDIR=/home/stewarz2/various_programs/uniref_db
DBFILE=uniref90.fasta

# Specify output details
OUTPREFIX=fortunella

# Specify tmp directory
TMPDIR=citrus_uniref_tmp

# Specify program parameters
CPUS=12
NUM_ITERS=1

#################################


# Step 1: Add MMSEQS to path
export PATH=${MMSEQDIR}/:$PATH

# Step 2: Run MMseqs2
python ${VARSCRIPTDIR}/run_mmseqs2.py -q ${QUERYDIR}/${QUERYFILE} -t ${DBDIR}/${DBFILE} -o ${OUTPREFIX} -m ${MMSEQDIR} -c ${CPUS} -n ${NUM_ITERS} --resume --blast_sort --tmp ${TMPDIR}
