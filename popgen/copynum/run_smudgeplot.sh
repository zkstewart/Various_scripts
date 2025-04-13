#!/bin/bash -l
#PBS -N smudgeplot
#PBS -l walltime=24:00:00
#PBS -l mem=80G
#PBS -l ncpus=8

cd $PBS_O_WORKDIR

conda activate smudgeplot

####

# Specify Various_scripts dir
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the FASTQ folder and suffix
FQDIR=/work/ePGL/resequencing/NGS_647_Pete_WGS/trimmed_reads
FQSUFFIX=P.fq.gz

# Specify the results directory
OUTDIR=smudgeplot_results

# Specify computational parameters
CPUS=8
MEM=60

####

# Run the pipeline
python ${VARSCRIPTDIR}/popgen/copynum/smudgeplot_pipeline.py \
    -i ${FQDIR} \
    -f fastq -k fastk \
    --fileSuffix ${FQSUFFIX} \
    -o ${OUTDIR} \
    --cpus ${CPUS} --mem ${MEM}
