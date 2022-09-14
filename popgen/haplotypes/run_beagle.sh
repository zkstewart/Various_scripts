#!/bin/bash -l
#PBS -N beagle
#PBS -l walltime=48:00:00
#PBS -l mem=80G
#PBS -l ncpus=8

cd $PBS_O_WORKDIR

module load java/1.8.0_231

#################################

# Specify the location of the Various_scripts directory
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the location of the Beagle JAR file
JARDIR=/home/stewarz2/various_programs/beagle
JARFILE=beagle.22Jul22.46e.jar

# Specify the location of the filtered VCF
VCFDIR=/home/stewarz2/plant_group/craig_analysis/almond/freebayes
VCFFILE=craig_almond.final.vcf

# Specify the location of the pop map file
POPSFILE=/home/stewarz2/plant_group/craig_analysis/almond/metadata/pop_map.txt

# Specify a prefix for output files
OUTPREFIX=craig_almond

# Specify computation resource parameters
CPUS=8

#################################

MINSNPS=5

# Step 1: Make VCF suitable for Beagle
PREFIX="${VCFFILE%.*}"
bcftools norm -D ${VCFDIR}/${VCFFILE} > ${PREFIX}.nodupe.vcf
python ${VARSCRIPTDIR}/popgen/VCF/filter_vcf.py -v ${PREFIX}.nodupe.vcf -p ${POPSFILE} -o ${PREFIX}.beagle.vcf --mpp 0 --min ${MINSNPS}

# Step 2: Run Beagle
java -Xmx64g -jar ${JARDIR}/${JARFILE} gt=${PREFIX}.beagle.vcf out=${OUTPREFIX} nthreads=${CPUS}
