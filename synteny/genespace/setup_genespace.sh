#!/bin/bash -l
#PBS -N setupGENESPACE
#PBS -l walltime=12:00:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify conda environments to use for different steps
PYENV=base
AGATENV=perl5

# Specify program/script locations
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts
GENSCRIPTDIR=/home/stewarz2/scripts/Genome_analysis_scripts
AGATDIR=/home/stewarz2/various_programs/AGAT

# Specify input file locations
P1=glauca
P2=australasica
P3=australis
P4=hindsii
P5=limon
P6=sinensis

F1DIR=/home/stewarz2/plant_group/glauca/genome
F2DIR=/home/stewarz2/plant_group/glauca/other_genomes
F3DIR=/home/stewarz2/plant_group/glauca/other_genomes
F4DIR=/home/stewarz2/plant_group/glauca/other_genomes
F5DIR=/home/stewarz2/plant_group/glauca/other_genomes
F6DIR=/home/stewarz2/plant_group/glauca/other_genomes

F1=C.glauca_collapsed_genome.fasta
F2=australasica_genome.fasta
F3=australis_genome.fasta
F4=hindsii_genome.fasta
F5=limon_genome.fasta
F6=sinensis_genome.fasta

A1=/home/stewarz2/plant_group/glauca/annotation/glauca_collapsed_braker.gff3
A2=/home/stewarz2/plant_group/glauca/other_genomes/australasica_annotation.gff3
A3=/home/stewarz2/plant_group/glauca/other_genomes/australis_annotation.gff3
A4=/home/stewarz2/plant_group/glauca/other_genomes/hindsii_annotation.gff3
A5=/home/stewarz2/plant_group/glauca/other_genomes/limon_annotation.gff3
A6=/home/stewarz2/plant_group/glauca/other_genomes/sinensis_annotation.gff3

# Specify files to flip
P4_TO_FLIP="chr1 chr3 chr6 chr7"
P5_TO_FLIP="GWHCBFQ00000001.1 GWHCBFQ00000005.1 GWHCBFQ00000011.1"
P6_TO_FLIP="chr2_v3.0 chr5_v3.0 chr7_v3.0 chr8_v3.0"

####

# STEP 0: Make sure genomes are multiline to avoid issues with AGAT
conda activate ${PYENV}
mkdir -p genomes

python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f single2multi -n 80 -i ${F1DIR}/${F1} -o genomes/${F1}
python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f single2multi -n 80 -i ${F2DIR}/${F2} -o genomes/${F2}
python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f single2multi -n 80 -i ${F3DIR}/${F3} -o genomes/${F3}
python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f single2multi -n 80 -i ${F4DIR}/${F4} -o genomes/${F4}
python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f single2multi -n 80 -i ${F5DIR}/${F5} -o genomes/${F5}
python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f single2multi -n 80 -i ${F6DIR}/${F6} -o genomes/${F6}

# STEP 1: Split apart GFF3 contigs that need reverse complementing ("flipping")
## >> 1.1: P4
rm -f ${P4}_toflip.txt
for VALUE in ${P4_TO_FLIP}; do echo ${VALUE} >> ${P4}_toflip.txt; done
python ${GENSCRIPTDIR}/gff3_entry_retrieve_remove.py -g ${A4} -t ${P4}_toflip.txt -m retrieve -b main -o ${P4}_toflip.gff3
python ${GENSCRIPTDIR}/gff3_entry_retrieve_remove.py -g ${A4} -t ${P4}_toflip.txt -m remove -b main -o ${P4}_noflip.gff3

## >> 1.2: P5
rm -f ${P5}_toflip.txt
for VALUE in ${P5_TO_FLIP}; do echo ${VALUE} >> ${P5}_toflip.txt; done
python ${GENSCRIPTDIR}/gff3_entry_retrieve_remove.py -g ${A5} -t ${P5}_toflip.txt -m retrieve -b main -o ${P5}_toflip.gff3
python ${GENSCRIPTDIR}/gff3_entry_retrieve_remove.py -g ${A5} -t ${P5}_toflip.txt -m remove -b main -o ${P5}_noflip.gff3

## >> 1.3: P6
rm -f ${P6}_toflip.txt
for VALUE in ${P6_TO_FLIP}; do echo ${VALUE} >> ${P6}_toflip.txt; done
python ${GENSCRIPTDIR}/gff3_entry_retrieve_remove.py -g ${A6} -t ${P6}_toflip.txt -m retrieve -b main -o ${P6}_toflip.gff3
python ${GENSCRIPTDIR}/gff3_entry_retrieve_remove.py -g ${A6} -t ${P6}_toflip.txt -m remove -b main -o ${P6}_noflip.gff3


# STEP 2: Flip GFF3s using AGAT
conda activate ${AGATENV}

## >> 2.1: P4
if [[ ! -f ${P4}_flipped.gff3 ]]; then
${AGATDIR}/bin/agat_sq_reverse_complement.pl --gff ${P4}_toflip.gff3 --fasta genomes/${F4} -o ${P4}_flipped.gff3;
fi;

## >> 2.2: P5
if [[ ! -f ${P5}_flipped.gff3 ]]; then
${AGATDIR}/bin/agat_sq_reverse_complement.pl --gff ${P5}_toflip.gff3 --fasta genomes/${F5} -o ${P5}_flipped.gff3;
fi;

## >> 2.3: P6
if [[ ! -f ${P6}_flipped.gff3 ]]; then
${AGATDIR}/bin/agat_sq_reverse_complement.pl --gff ${P6}_toflip.gff3 --fasta genomes/${F6} -o ${P6}_flipped.gff3;
fi;


# STEP 3: Join GFF3s back together
mkdir -p setup_gff3s

## >> 3.1: P4
cat ${P4}_noflip.gff3 ${P4}_flipped.gff3 > setup_gff3s/${P4}_prepared.gff3

## >> 3.2: P5
cat ${P5}_noflip.gff3 ${P5}_flipped.gff3 > setup_gff3s/${P5}_prepared.gff3

## >> 3.3: P6
cat ${P6}_noflip.gff3 ${P6}_flipped.gff3 > setup_gff3s/${P6}_prepared.gff3


# STEP 4: Generate BED format from GFF3
conda activate ${PYENV}
mkdir -p bed

## >> 4.1: Default
python ${VARSCRIPTDIR}/GFF3/gff3_id_mapper.py -forEach gene -map contig -to start end ID -g ${A1} -o bed/${P1}.bed --noHeader --tolerant
python ${VARSCRIPTDIR}/GFF3/gff3_id_mapper.py -forEach gene -map contig -to start end ID -g ${A2} -o bed/${P2}.bed --noHeader --tolerant
python ${VARSCRIPTDIR}/GFF3/gff3_id_mapper.py -forEach gene -map contig -to start end ID -g ${A3} -o bed/${P3}.bed --noHeader --tolerant

## >> 4.2: Flipped
python ${VARSCRIPTDIR}/GFF3/gff3_id_mapper.py -forEach gene -map contig -to start end ID -g setup_gff3s/${P4}_prepared.gff3 -o bed/${P4}.bed --noHeader --tolerant
python ${VARSCRIPTDIR}/GFF3/gff3_id_mapper.py -forEach gene -map contig -to start end ID -g setup_gff3s/${P5}_prepared.gff3 -o bed/${P5}.bed --noHeader --tolerant
python ${VARSCRIPTDIR}/GFF3/gff3_id_mapper.py -forEach gene -map contig -to start end ID -g setup_gff3s/${P6}_prepared.gff3 -o bed/${P6}.bed --noHeader --tolerant


# STEP 5: Generate peptide FASTA files
## Note that we don't use the flipped GFF3s since they aren't guaranteed to have splice frames set right by AGAT
if [[ ! -f ${P1}.aa ]]; then
python ${GENSCRIPTDIR}/gff3_to_fasta_v2.py -i ${F1DIR}/${F1} -g ${A1} -l main -s cds -o ${P1} --relaxed --gene_ids;
fi;

if [[ ! -f ${P2}.aa ]]; then
python ${GENSCRIPTDIR}/gff3_to_fasta_v2.py -i ${F2DIR}/${F2} -g ${A2} -l main -s cds -o ${P2} --relaxed --gene_ids;
fi;

if [[ ! -f ${P3}.aa ]]; then
python ${GENSCRIPTDIR}/gff3_to_fasta_v2.py -i ${F3DIR}/${F3} -g ${A3} -l main -s cds -o ${P3} --relaxed --gene_ids;
fi;

if [[ ! -f ${P4}.aa ]]; then
python ${GENSCRIPTDIR}/gff3_to_fasta_v2.py -i ${F4DIR}/${F4} -g ${A4} -l main -s cds -o ${P4} --relaxed --gene_ids;
fi;

if [[ ! -f ${P5}.aa ]]; then
python ${GENSCRIPTDIR}/gff3_to_fasta_v2.py -i ${F5DIR}/${F5} -g ${A5} -l main -s cds -o ${P5} --relaxed --gene_ids;
fi;

if [[ ! -f ${P6}.aa ]]; then
python ${GENSCRIPTDIR}/gff3_to_fasta_v2.py -i ${F6DIR}/${F6} -g ${A6} -l main -s cds -o ${P6} --relaxed --gene_ids;
fi;

# STEP 6: Symlink just the peptide files to their respective folder with GENESPACE-expected file names
mkdir -p peptide
ln -s ../${P1}.aa peptide/${P1}.fa
ln -s ../${P2}.aa peptide/${P2}.fa
ln -s ../${P3}.aa peptide/${P3}.fa
ln -s ../${P4}.aa peptide/${P4}.fa
ln -s ../${P5}.aa peptide/${P5}.fa
ln -s ../${P6}.aa peptide/${P6}.fa

# STEP 7: Now run R!
