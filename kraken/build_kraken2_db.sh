#!/bin/bash -l
#PBS -N k2DB
#PBS -l walltime=48:00:00
#PBS -l mem=300G
#PBS -l ncpus=12

cd $PBS_O_WORKDIR

####

# Specify the kraken conda environment
CONDAENV=kraken2

# Specify the annotarium GitHub repository location
ANNOTARIUM=/home/stewarz2/scripts/annotarium

# Specify the database directory and temp files location
TMPDIR=/work/ePGL/databases/kraken2/tmp
DBDIR=/work/ePGL/databases/kraken2/citrus_plants_contaminants

# Specify the location of additional genomes to add into the database
EXTRA1=/work/ePGL/genomes/citrus/aurantiifolia/mxlime_usda_v1/GCA_052148995.1_Mxlime_USDA_v1_genomic.fna
EXTRA2=/work/ePGL/genomes/citrus/australasica/rainbow_v1/rainbow_v1_hap1.fasta
EXTRA3=/work/ePGL/genomes/citrus/australis/uq_v1_nakandala/australis.fasta
EXTRA4=/work/ePGL/genomes/citrus/clementina/GCA_037179405.1/GCA_037179405.1_UCR_Fairchild_v1.0_genomic.fna
EXTRA5=/work/ePGL/genomes/citrus/hindsii/gj_v2/GJ.v2.0.genome.fa
EXTRA6=/work/ePGL/genomes/citrus/kinokuni/kazusa_kishu_mikan/kinokuni_kz_hap1.fasta
EXTRA7=/work/ePGL/genomes/citrus/maxima/guangdong_zm/maxima_xm_hap1.fasta
EXTRA8=/work/ePGL/genomes/citrus/murcott/henry_upuli_2025/murcott_hap1.fasta
EXTRA9=/work/ePGL/genomes/citrus/reticulata/ucsk1/reticulata_ucsk1_hap1.fasta

# Specify the taxid of each additional genome
TAXID1=159033
TAXID2=416196
TAXID3=341934
TAXID4=85681
TAXID5=159041
TAXID6=408488
TAXID7=37334
TAXID8=197943
TAXID9=85571

# Specify number of threads to use
CPUS=12

####

# Modify the sequence headers of additional genomes for database inclusion
python ${ANNOTARIUM}/annotarium.py fasta rename -i ${EXTRA1} -o ${TMPDIR}/${TAXID1}.fasta --format "{seqid}|kraken:taxid|${TAXID1}"
python ${ANNOTARIUM}/annotarium.py fasta rename -i ${EXTRA2} -o ${TMPDIR}/${TAXID2}.fasta --format "{seqid}|kraken:taxid|${TAXID2}"
python ${ANNOTARIUM}/annotarium.py fasta rename -i ${EXTRA3} -o ${TMPDIR}/${TAXID3}.fasta --format "{seqid}|kraken:taxid|${TAXID3}"
python ${ANNOTARIUM}/annotarium.py fasta rename -i ${EXTRA4} -o ${TMPDIR}/${TAXID4}.fasta --format "{seqid}|kraken:taxid|${TAXID4}"
python ${ANNOTARIUM}/annotarium.py fasta rename -i ${EXTRA5} -o ${TMPDIR}/${TAXID5}.fasta --format "{seqid}|kraken:taxid|${TAXID5}"
python ${ANNOTARIUM}/annotarium.py fasta rename -i ${EXTRA6} -o ${TMPDIR}/${TAXID6}.fasta --format "{seqid}|kraken:taxid|${TAXID6}"
python ${ANNOTARIUM}/annotarium.py fasta rename -i ${EXTRA7} -o ${TMPDIR}/${TAXID7}.fasta --format "{seqid}|kraken:taxid|${TAXID7}"
python ${ANNOTARIUM}/annotarium.py fasta rename -i ${EXTRA8} -o ${TMPDIR}/${TAXID8}.fasta --format "{seqid}|kraken:taxid|${TAXID8}"
python ${ANNOTARIUM}/annotarium.py fasta rename -i ${EXTRA9} -o ${TMPDIR}/${TAXID9}.fasta --format "{seqid}|kraken:taxid|${TAXID9}"

# Build the database with inclusion of relevant online databases
conda activate ${CONDAENV}

kraken2-build --threads ${CPUS} --use-ftp --standard --db ${DBDIR} # Refseq archaea, bacteria, viral, plasmid, human, UniVec_Core
kraken2-build --threads ${CPUS} --use-ftp --download-library fungi --db ${DBDIR}
kraken2-build --threads ${CPUS} --use-ftp --download-library plant --db ${DBDIR}

# Add in extra local genomes for databasing
kraken2-build --threads ${CPUS} --add-to-library ${TMPDIR}/${TAXID1}.fasta --db ${DBDIR}
kraken2-build --threads ${CPUS} --add-to-library ${TMPDIR}/${TAXID2}.fasta --db ${DBDIR}
kraken2-build --threads ${CPUS} --add-to-library ${TMPDIR}/${TAXID3}.fasta --db ${DBDIR}
kraken2-build --threads ${CPUS} --add-to-library ${TMPDIR}/${TAXID4}.fasta --db ${DBDIR}
kraken2-build --threads ${CPUS} --add-to-library ${TMPDIR}/${TAXID5}.fasta --db ${DBDIR}
kraken2-build --threads ${CPUS} --add-to-library ${TMPDIR}/${TAXID6}.fasta --db ${DBDIR}
kraken2-build --threads ${CPUS} --add-to-library ${TMPDIR}/${TAXID7}.fasta --db ${DBDIR}
kraken2-build --threads ${CPUS} --add-to-library ${TMPDIR}/${TAXID8}.fasta --db ${DBDIR}
kraken2-build --threads ${CPUS} --add-to-library ${TMPDIR}/${TAXID9}.fasta --db ${DBDIR}

# Build the db
kraken2-build --threads ${CPUS} --build --db ${DBDIR}
