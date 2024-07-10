#!/bin/bash -l
#PBS -N cafepipe
#PBS -l ncpus=8
#PBS -l walltime=48:00:00
#PBS -l mem=20G

cd $PBS_O_WORKDIR

conda activate cafe

####

# Specify Various_scripts location
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify input file locations
COUNTS=/home/stewarz2/plant_group/glauca/cafe/selfliftoff/orthofinder/analysis/OrthoFinder/Results_Jul09/Orthogroups/Orthogroups.GeneCount.tsv
TREE=/home/stewarz2/plant_group/glauca/cafe/selfliftoff/tree/glauca_sc_gene_msa.fasta.timetree.nwk

# Specify new sequence IDs (set during singlecopy MSA build)
NEWIDS="australasica australis glauca hindsii limon sinensis"

# Specify computation parameters
CPUS=8
PERCENTILE=99.95 # 99.95 filters the top 0.5 percent of families

####

# STEP 1: Modify GeneCount table for CAFE5 use
python ${VARSCRIPTDIR}/Orthofinder/cafe/orthogroup_genecount_to_cafe.py \
    -i ${COUNTS} \
    -o Orthogroups.GeneCount.CAFE5.tsv \
    --new_ids ${NEWIDS}

# STEP 2: Filter the file for CAFE input
python ${VARSCRIPTDIR}/CAFE/filter_families_on_difference.py \
    -i Orthogroups.GeneCount.CAFE5.tsv \
    -o Orthogroups.GeneCount.CAFE5.filtered.tsv \
    --percentile ${PERCENTILE}

# STEP 3: Run CAFE
for gammarate in $(seq 1 5); do
    cafe5 -i Orthogroups.GeneCount.CAFE5.filtered.tsv -t ${TREE} -c ${CPUS} -k ${gammarate} -o gamma${gammarate};
done;


########################################
# DOWNSTREAM PROCESSING OF BEST RESULT #
########################################


# STEP 4: Get multiple runs of optimal gamma value
gammarate=2
for i in $(seq 1 5); do
    cafe5 -i Orthogroups.GeneCount.CAFE5.filtered.tsv -t ${TREE} -c ${CPUS} -k ${gammarate} -o gamma${gammarate}_run${i};
done;

# STEP 5: Plot the chosen run with the smallest likelihood value
chosen=gamma2
cafeplotter -i ${chosen} -o ${chosen}_plots --ignore_branch_length --format pdf;
