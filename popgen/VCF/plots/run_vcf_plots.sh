#!/bin/bash -l
#PBS -N vcfplot
#PBS -l walltime=32:00:00
#PBS -l mem=70G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify the location of the variantopia github repository
VARIANTOPIA=/home/stewarz2/scripts/variantopia

# Specify the genome FASTA file location
GENOME=/home/stewarz2/flgenomics/genomes/filename.fasta

# Specify the VCF file location
VCF=/home/stewarz2/flgenomics/commercial/variant_calls/merged.filtered.vcf.gz

# Specify the prefix for output files
PREFIX=flgenomics_commercial

# Specify behavioural parameters
OUTFORMAT=png # png OR svg OR pdf
COLOUR=viridis # viridis OR Greys OR GnBu OR RdBu; see https://matplotlib.org/stable/users/explain/colors/colormaps.html#sequential
WINDOWSIZE=100000 # length of genome to summarise statistics within
WIDTH=10 # plot dimension
HEIGHT=6 # plot dimension

####

for STATISTIC in "snpnumber" "maf" "callrate" "het"; do
    python ${VARIANTOPIA}/variantopia.py \
	    vcf plot -i ${VCF} \
                 -o ${PREFIX}.${STATISTIC}.${OUTFORMAT} \
                 -s ${STATISTIC} \
                 -f chromosomes \
                 -w ${WINDOWSIZE} \
                 --genome ${GENOME} \
                 --colour ${COLOUR} \
                 --width ${WIDTH} --height ${HEIGHT};
done
