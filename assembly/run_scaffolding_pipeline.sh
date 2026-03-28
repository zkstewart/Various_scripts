#!/bin/bash -l
#PBS -N scaffolding
#PBS -l walltime=12:00:00
#PBS -l mem=70G
#PBS -l ncpus=8

cd $PBS_O_WORKDIR

#################################

# Specify the location of the Various_scripts repository
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the location of the diploid reference genome
REFDIR=/work/ePGL/genomes/citrus/maxima/guangdong_zm
REFHAP1=maxima_xm_hap1.rename.fasta
REFHAP2=maxima_xm_hap2.rename.fasta

# Specify how many CPUs to use
CPUS=8

# Specify behavioural params
OUTDIR=output
PRESET=asm20

#################################

# STEP 0: Specify the location of the input files
ASSEMBLYDIR=$(realpath ../assembly)
HAP1=$(ls ${ASSEMBLYDIR}/*.hap1.p_ctg.fasta)
HAP2=$(ls ${ASSEMBLYDIR}/*.hap2.p_ctg.fasta)

# STEP 1: Run the pipeline
python ${VARSCRIPTDIR}/assembly/haplotype_scaffolding_pipeline.py -i ${HAP1} ${HAP2} \
    -r ${REFDIR}/${REFHAP1} ${REFDIR}/${REFHAP2} \
    -o ${OUTDIR} --threads ${CPUS} --preset ${PRESET}
