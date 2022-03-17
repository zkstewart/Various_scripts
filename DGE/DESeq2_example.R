# DESeq2_example.R
# Example script for running DESeq2 in various modes

# Helpful websites
## Overall DESeq2 guide
# http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments

## Time-course analyses help
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
# https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

# Install packages
## The below libraries need to be installed from Bioconductor / CRAN
## Uncomment (remove the leading "# " from) the below lines if you need to do this

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("DEGreport")
# BiocManager::install("goseq")
# 
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("pheatmap")
# install.packages("glmpca")
# install.packages("RColorBrewer")

# Load packages
library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(glmpca)
library(RColorBrewer)
# library(goseq) # Not detailing goseq here
# library(DEGreport) # Not detailing clustering here

# Setup directory locations
## To re-use this code most effectively, it helps to have all your file name
## variables declared up here. That way you know exactly where to change things
## when you run this on another dataset and you don't need to mess with the code
## that's doing things.
## NOTE: Some of the file names still need to be set down below, but you shouldn't
## need to change code lines much. Just look for capitalised variables.

# >> DGE-related variables
## Note that directory slashes need to go forwards ("/") even on a Windows computer
COUNTS_FILE = "F:/plant_group/plant_rnaseq/DEW_counts/alm/alm_dew_counts.fix.tsv"
COLDATA_FILE = "F:/plant_group/plant_rnaseq/sample_metadata/alm/alm_coldata.txt"



# >> goseq-related variables
## If you're not planning to run goseq, this is irrelevant
## This script currently doesn't perform goseq, so you can ignore it!
#GOANNOT_FILE = "F:/plant_group/plant_rnaseq/annotations/alm/DEW_annotations/alm.GO.fix.mapping"
#LENGTHS_FILE = "F:/plant_group/plant_rnaseq/annotations/alm/DEW_annotations/Prudul26A.gene.fix.lengths"



# >> Working directory variable
WORK_DIR = "F:/testing" # Set this to be the location of this script file, where outputs will be generated
setwd(WORK_DIR)



# Create output directories
## You shouldn't need to change this - goseq results will be stored
## in this folder.
CLUSTER_DIR = paste0(WORK_DIR, "/clustering")
GOSEQ_DIR = paste0(WORK_DIR, "/goseq")

dir.create(CLUSTER_DIR, showWarnings = FALSE)
dir.create(GOSEQ_DIR, showWarnings = FALSE)



# Specify significance threshold
## Sometimes you want your DGE and goseq to run at different thresholds.
## Set these to values that give you a manageable amount of results (keep it <=0.05 though)
SIGNIFICANCE_CUTOFF = 0.001
GOSEQ_SIGNIFICANCE_THRESHOLD = 0.05



# Load in data files for DGE analysis
# >> Counts file
## If you've got a table like the below...
## +-----------------+------+------+-----+
## |     gene_id     | alm1 | alm2 | ... |
## +-----------------+------+------+-----+
## | Prudul26A002130 |  166 |  128 | ... |
## | Prudul26A026237 |  260 |  118 | ... |
## +-----------------+------+------+-----+
## ... then you should set row.names to be "gene_id"
counts.table = read.table(file=COUNTS_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = "gene_id")



# >> Coldata file
## If you've got a table like the below...
## +------------+--------+------+-----------+
## | SampleName | tissue | time | replicate |
## +------------+--------+------+-----------+
## | alm1       | L      |    1 |         1 |
## | alm2       | L      |    1 |         2 |
## | alm3       | L      |    2 |         1 |
## | alm4       | L      |    2 |         2 |
## | alm5       | B      |    1 |         1 |
## | alm6       | B      |    1 |         2 |
## | alm7       | B      |    2 |         1 |
## | alm8       | B      |    2 |         2 |
## +------------+--------+------+-----------+
## ... then you should set row.names to be "SampleName"
coldata.table = read.table(file=COLDATA_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = "SampleName")



# Specify as factor any relevant numeric variables that will be part of the experimental design
## Any number that will be used as a variable in the DGE analysis needs to be set to be a factor.
## Even time variables in a time course analysis should be made a factor for various reasons that
## go beyond my understanding. Just make them a factor. Variables you're NOT using in the analysis
## but are numbers e.g., replicates, can be left alone.
coldata.table$time <- factor(coldata.table$time)




###########################################################################
##                                                                       ##
##                        DATA EXPLORATION                               ##
##                                                                       ##
###########################################################################



## In this section, we're not running a DGE. But we are using DESeq2's counts
## values to do some overall exploration of the data.

# Create a general purpose DESeq2 object for exploration
## No changes will be needed here.
dds <- DESeqDataSetFromMatrix(countData = counts.table, colData = coldata.table, design = ~ 1)
dds <- estimateSizeFactors(dds)

# Output normalised counts
## Set this file name to something appropriate. This is the table you should
## use when creating expression plots.
NORM_COUNTS_FILE = "my.normalised.counts"
write.table(as.data.frame(counts(dds, normalized=TRUE)), file=NORM_COUNTS_FILE, sep="\t", quote=FALSE, col.names=NA)

# Drop low count genes
## No changes needed here.
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 2
dds <- dds[filter,]


#### Distances heatmap

# Calculate distances
## vst() performs a variance stabilising transformation. I don't get the maths, but it's necessary
## when considering the notion of "distance" between samples.
vsd <- vst(dds, blind = FALSE) 
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )

# Create heatmap
## The below rownames() line should contain a combination of variables to uniquely identify your sample.
## This creates labels for the heatmap. Set "tissue" to be "treatment" for example if that's relevant.
rownames(sampleDistMatrix) <- paste( vsd$tissue, vsd$replicate, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors) # Need to save this manually in R studio with export

#### PCA
## In this PCA, we're actually performing a GLM-PCA, which is a fancy form of PCA
## that is probably better at handling data with various biases. It's just more
## robust and I opt for it when I can.

# Run the main GLM-PCA computation
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors

# Associate any variables you want plotted to the table
gpca.dat$tissue <- dds$tissue
gpca.dat$time <- dds$time
gpca.dat$replicate <- dds$replicate

# Create the plot
## You might want to play with the color and shape settings, based on the variables
## you associated above. Pick two variables that explain your data well. If you need to
## show more variables somehow, you might need to combine variables or facet your plot somehow.
PCA_FILE_NAME="alm_glmpca.png"

png(filename=PCA_FILE_NAME, width = 480, height = 480) # play with width and height values to make it look good
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = time, shape = tissue)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA") # play with the Title, and any other ggplot stuff
dev.off()




###########################################################################
##                                                                       ##
##                  DIFFERENTIAL GENE EXPRESSION                         ##
##                                                                       ##
###########################################################################

## This section will show two ways to run DGE with DESeq2 i.e., pairwise
## comparison and time-course analysis.


####### METHOD 1: Pairwise comparison

# Create a DESeq2 object for pairwise comparisons
## It's important to consider your design here, but it should be easy to set it.
## Simply specify all the variables of relevance to your experiment with "+" in
## between them. Don't include irrelevant variables like replicate. If this isn't
## giving you what you want, you might need to combine variables in your colData file
## so you have a single "group" variable or something of the sort.
dds <- DESeqDataSetFromMatrix(countData = counts.table, colData = coldata.table, design = ~ tissue + time)
dds <- estimateSizeFactors(dds)

# Drop low count genes
## No changes needed here.
## For your interest, doing this step can be important to prevent
## low quality results from being flagged as DE just because one or
## two libraries showed minimal expression when all others did not.
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 2
dds <- dds[filter,]

# Run DESeq2
## Nothing should need changing here.
dds <- DESeq(dds)
res <- results(dds, alpha=SIGNIFICANCE_CUTOFF)
res <- res[order(res$padj),]
res <- na.omit(res)

# Extract pairwise comparison results
## Nothing should need changing here.
comparisons = resultsNames(dds)[2:length(resultsNames(dds))]
for (c in comparisons)
{
  # Extract relevant details from comparison
  cSplit = strsplit(c, "_")
  contrast = cSplit[[1]][1]
  group1 = cSplit[[1]][2]
  group2 = cSplit[[1]][4]

  # Obtain pairwise comparison of interest
  res.pw <- results(dds, alpha=SIGNIFICANCE_CUTOFF, contrast=c(contrast, group1, group2))
  
  # Clean up
  res.pw <- res.pw[order(res.pw$padj),]
  res.pw <- na.omit(res.pw)
  
  # Get DE genes
  res.DE.pw <- res.pw[res.pw$padj < SIGNIFICANCE_CUTOFF,]
  
  # Output table
  fileName = paste0(c, "_deseq2_results.tsv")
  write.table(as.data.frame(res.DE.pw), file=fileName, sep="\t", quote=FALSE, col.names=NA)
}
## After this, you should have a bunch of files like "time_6_vs_1_deseq2_results.tsv"
## Included are only results that pass our SIGNIFICANCE_CUTOFF value.



####### METHOD 2: Time-course analysis

# Create a DESeq2 object for time-course analysis
## It's important to consider your design here, but it should be easy to set it once again.
## This design should only include your time variable, whatever it's called. If you have multiple
## tissues/treatments you want to run as time courses separately, you can do it in R, but I'd 
## recommend just subsetting your countData and colData files to only include the samples you
## want to analyse.
dds.tc <- DESeqDataSetFromMatrix(countData = counts.table, colData = coldata.table, design = ~ time)
dds.tc <- estimateSizeFactors(dds.tc)

# Drop low count genes
## Same as before
nc.tc <- counts(dds.tc, normalized=TRUE)
filter.tc <- rowSums(nc.tc >= 10) >= 2
dds.tc <- dds.tc[filter.tc,]

# Run DESeq2
## Nothing should need changing here.
## We use a likelihood ratio test (LRT) with a reduced model.
## Long story short, if you've subsetted data so you only have your treatment of
## interest across its time points, a ~1 model will work. If you have other
## "confounding" variables you kind of want to "control" for, you'd specify that
## as your reduced model e.g., using "~tissue" would essentially tell DESeq2
## to find genes that change over time controlling for any differences we see between
## the tissues. But, usually, you'd rather have the time course for the tissues/treatments
## separately.
dds.tc <- DESeq(dds.tc, test = "LRT", reduced = ~1)
res.tc <- results(dds.tc, alpha=SIGNIFICANCE_CUTOFF)
res.tc <- res.tc[order(res.tc$padj),]
res.tc <- na.omit(res.tc)

# Get DE genes
res.DE.tc <- res.tc[res.tc$padj < SIGNIFICANCE_CUTOFF,]

# Output table
## Choose an appropriate file name to write to
TC_RESULTS_NAME = paste0("leafTissue_deseq2_TC_results.tsv")
write.table(as.data.frame(res.DE.tc), file=TC_RESULTS_NAME, sep="\t", quote=FALSE, col.names=NA)


###########################################################################
##                                                                       ##
##                       DGE RESULTS PLOTTING                            ##
##                                                                       ##
###########################################################################

## Once you've got your results, there's a few things you might want to do.
## This script will only show two basic plots - boxplots and a DE gene heatmap.


####### METHOD 1: Boxplots

# Create the directory where plots will be saved
PLOTS_DIR = "expression_boxplots"
dir.create(PLOTS_DIR)

# Get the list of DE genes
## Here, I'm using "res.DE.tc" which we calculated as part of the time-course
## analysis. But you can put "res.DE.pw" or any DESeq2 table after you've
## done the SIGNIFICANT_CUTOFF filtering of results.
gene.ids = rownames(res.DE.tc)

# Get a table of normalised counts
## When plotting expression, a fair comparison of different RNAseq libraries
## can only occur after data normalisation. It helps to deal with biases
## in sequencing depth as well as the impact that gene upregulation has
## on all other genes (e.g., if gene A goes up in expression, we expect to
## sequence less of every other gene; this isn't biologically real, but it
## is reflected in RNA sequencing).
normCountsTable = as.data.frame(counts(dds.tc, normalized=TRUE)) # Use the "dds" you performed the analysis with

# Generate plots for all DE genes
## Make the plots look different by changing the ggplot line.
## Otherwise, make sure to specify an appropriate row$group value from
## your colData file as the comment below asks for. The variable doesn't
## have to be called "group" either. It's up to you to make this look good
## and be appropriate for your analysis.
for (g.id in gene.ids) {
  
  row = normCountsTable[rownames(normCountsTable) == g.id,]
  row = as.data.frame(t(row))
  colnames(row) = c("expression")
  
  # Change the "group" value to be a variable of interest from your colData file!
  row$group = coldata.table$time
  
  p = ggplot(row, aes(x = group, y=expression, fill=group)) +
    geom_boxplot() +
    theme_bw()
  
  ggsave(paste0(PLOTS_DIR, "/", g.id, ".png"))
}

####### METHOD 2: Heatmap

# Function for computing Z-score
## Don't touch this.
## Function credit to https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Create the directory where plots will be saved
PLOTS_2_DIR = "expression_heatmaps"
dir.create(PLOTS_2_DIR)

# Get the list of DE genes
## Refer to method 1 comments for this process.
gene.ids = rownames(res.DE.tc)

# Get a table of normalised counts
## Again, refer to method 1. This works a bit differently since we're
## subsetting it to only have the counts for our DE genes.
normCountsTable = as.data.frame(counts(dds.tc, normalized=TRUE)) # Use the "dds" you performed the analysis with
normCountsTable.DE = normCountsTable[rownames(normCountsTable) %in% gene.ids,]

# Get our plot labels
## This works the same as in the distance heatmap detailed in the exploration section.
## Make sure your "dds" is the one you want to use, and choose appropriate variables
## that uniquely identify each sample.
colnames(normCountsTable.DE) <- paste( dds.tc$tissue, dds.tc$time, sep = " - " )
normCountsTable.DE.Znorm <- t(apply(normCountsTable.DE, 1, cal_z_score))
pheatmap(normCountsTable.DE.Znorm) # Need to save this manually in R studio with export




###########################################################################
##                                                                       ##
##                           FINAL NOTES                                 ##
##                                                                       ##
###########################################################################

# There's things I haven't included in this R script, like how to perform
# clustering. It's not difficult, but it can take a bit more configuring to
# get it right and you might see errors if there's not enough DE genes coming
# out of your analysis.

# If you want further help, email me at zkstewart1@gmail.com.


