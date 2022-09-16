#!/bin/bash -l
#PBS -N beagleH
#PBS -l walltime=04:00:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

#################################

# Specify the location of the Various_scripts directory
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the location of the Beagle VCF
VCFDIR=/home/stewarz2/plant_group/craig_analysis/citrus/beagle
VCF=craig_citrus.vcf

# Specify the location of the genome FASTA
FASTADIR=/home/stewarz2/plant_group/craig_analysis/citrus/genome
FASTA=clementina.fasta

# Specify the location of the genome annotation GFF3
GFF3DIR=/home/stewarz2/plant_group/craig_analysis/citrus/annotation
GFF3=Cclementina_182_v1.0.gene_exons.gff3

# Specify the output directory
OUTDIR=haplotype_results

# Specify the minimum number of SNPs in a gene to be considered a haplotype
MINSNPS=3

#################################


# RUN
python ${VARSCRIPTDIR}/popgen/haplotypes/beagle_gene_haplotypes.py -v ${VCFDIR}/${VCF} \
	-g ${GFF3DIR}/${GFF3} \
	-f ${FASTADIR}/${FASTA} \
	-o ${OUTDIR} \
	--minSnps ${MINSNPS}
