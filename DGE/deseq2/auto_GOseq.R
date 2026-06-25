# BINge_GOseq.R

library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(goseq)
library(RColorBrewer)
library(comprehenr)
library(vegan)
library(grid)
library(gridExtra)
library(readr)
library(data.table)


###########################################################################
##                                                                       ##
##                     FUNCTION DECLARATIONS                             ##
##                                                                       ##
###########################################################################

parse_deseq2_degs <- function(fileName) {
  dge.table <- read.table(fileName, as.is=T, header=T, sep='\t')
  dge.ids <- dge.table$X
  return (dge.ids)
}

parse_wgcna_colours <- function(fileName) {
  colours.table <- read.table(fileName, as.is=T, header=T, sep='\t')
  return (colours.table)
}


###########################################################################
##                                                                       ##
##                         ANALYSIS SETUP                                ##
##                                                                       ##
###########################################################################

# Setup working directory variable
WORK_DIR = "F:/citrus/ted/DGE"
setwd(WORK_DIR)

# Specify significance threshold
GOSEQ_SIGNIFICANCE_THRESHOLD = 0.05

# Specify goseq-related variables
GOANNOT_FILE = "F:/citrus/ted/annotation/binge.gos.tsv"
LENGTHS_FILE = "F:/citrus/ted/binge/BINge_clustering_representatives.cds.lengths.tsv"

# Specify the DESeq2 output file location
DGE_DIR = "F:/citrus/ted/DGE/deseq2"

# Specify the WGCNA output file locations
WGCNA_BUD = "F:/citrus/ted/DGE/WGCNA_results/gene_colours.bud.tsv"
WGCNA_LEAF = "F:/citrus/ted/DGE/WGCNA_results/gene_colours.leaf.tsv"

###########################################################################
##                                                                       ##
##                             GOSEQ SETUP                               ##
##                                                                       ##
###########################################################################

# Generate GO mapping structure
goannot=read.table(GOANNOT_FILE, as.is=T, header=F, sep='\t')
goannot[goannot == 0] <- NA
splitgoannot = strsplit(goannot[,2], split='; ')
name.list = as.vector(goannot[,1])
names(splitgoannot) = name.list
rm(goannot) # tidy up

# Load in gene lengths as a vector
genelengths=read.table(LENGTHS_FILE, as.is=T, header=F)
prevector = genelengths[2]
rownames(prevector) = genelengths[,1] 
length.vector = as.integer(prevector[,1])
names(length.vector) = rownames(prevector)
rm(prevector) # tidy up
rm(genelengths) # tidy up

# Subset vectors to those shared in common
length.vector = length.vector[names(length.vector) %in% names(splitgoannot)]
splitgoannot = splitgoannot[names(splitgoannot) %in% names(length.vector)]

# Ensure consistency of length and GO data
length.names <- names(length.vector)
go.names <- names(splitgoannot)
if (! all(length.names %in% go.names))
{
  stop("Inconsistency 1 in length and GO files; fix this!")
}
if (! all(length.names == go.names))
{
  stop("Inconsistency 2 in metadata and counts; fix this!")
}

gene.names <- go.names # length and go are the same so doesn't matter which we pick

###########################################################################
##                                                                       ##
##                             GOSEQ: DEGs                               ##
##                                                                       ##
###########################################################################

# Create output directory
GOSEQ_DIR = paste0(WORK_DIR, "/goseq")
dir.create(GOSEQ_DIR, showWarnings = FALSE)

# Locate all DESeq2 result files
deseq2.files <- list.files(path = DGE_DIR, pattern = "\\.tsv$", full.names = TRUE)

# Run DESeq2 for each file
for (fileName in deseq2.files)
{
  # Determine output file name and skip if already generated
  output.file.name <- file.path(GOSEQ_DIR, basename(fileName))
  if (file.exists(output.file.name))
  {
    next
  }
  
  # Obtain the DE gene IDs
  deg.ids <- parse_deseq2_degs(fileName)
  
  # Create a vector identifying which genes are identified as DE in this analysis
  DE.vector <- gene.names %in% deg.ids
  DE.vector = DE.vector * 1 # converts TRUE to 1 and FALSE to 0
  DE.vector = as.integer(DE.vector)
  names(DE.vector) = gene.names
  
  # Run PWF for length bias
  pwf <- nullp(DEgenes=DE.vector, bias.data=length.vector)
  
  # Run goseq
  go_output <- goseq(pwf, gene2cat = splitgoannot)
  
  # Limit output to significant results
  go_output.SIG = go_output[go_output$over_represented_pvalue <= GOSEQ_SIGNIFICANCE_THRESHOLD | go_output$under_represented_pvalue <= GOSEQ_SIGNIFICANCE_THRESHOLD,]
  
  # Produce output table
  write.table(go_output.SIG, file=output.file.name, sep="\t", row.names = FALSE, quote=FALSE)
}

###########################################################################
##                                                                       ##
##                            GOSEQ: WGCNA                               ##
##                                                                       ##
###########################################################################

colours.bud = parse_wgcna_colours(WGCNA_BUD)
colours.leaf = parse_wgcna_colours(WGCNA_LEAF)

for (colour in unique(colours.bud$module))
{
  # Determine output file name and skip if already generated
  output.file.name <- file.path(GOSEQ_DIR, paste0(colour, ".bud.tsv"))
  if (file.exists(output.file.name))
  {
    next
  }
  
  # Obtain the module gene IDs
  colour.ids <- colours.bud[colours.bud$module == colour,]$gene_id
  
  # Create a vector identifying which genes are identified as DE in this analysis
  DE.vector <- gene.names %in% colour.ids
  DE.vector = DE.vector * 1 # converts TRUE to 1 and FALSE to 0
  DE.vector = as.integer(DE.vector)
  names(DE.vector) = gene.names
  
  # Run PWF for length bias
  pwf <- nullp(DEgenes=DE.vector, bias.data=length.vector)
  
  # Run goseq
  go_output <- goseq(pwf, gene2cat = splitgoannot)
  
  # Limit output to significant results
  go_output.SIG = go_output[go_output$over_represented_pvalue <= GOSEQ_SIGNIFICANCE_THRESHOLD | go_output$under_represented_pvalue <= GOSEQ_SIGNIFICANCE_THRESHOLD,]
  
  # Produce output table
  write.table(go_output.SIG, file=output.file.name, sep="\t", row.names = FALSE, quote=FALSE)
}

for (colour in unique(colours.leaf$module))
{
  # Determine output file name and skip if already generated
  output.file.name <- file.path(GOSEQ_DIR, paste0(colour, ".leaf.tsv"))
  if (file.exists(output.file.name))
  {
    next
  }
  
  # Obtain the module gene IDs
  colour.ids <- colours.leaf[colours.leaf$module == colour,]$gene_id
  
  # Create a vector identifying which genes are identified as DE in this analysis
  DE.vector <- gene.names %in% colour.ids
  DE.vector = DE.vector * 1 # converts TRUE to 1 and FALSE to 0
  DE.vector = as.integer(DE.vector)
  names(DE.vector) = gene.names
  
  # Run PWF for length bias
  pwf <- nullp(DEgenes=DE.vector, bias.data=length.vector)
  
  # Run goseq
  go_output <- goseq(pwf, gene2cat = splitgoannot)
  
  # Limit output to significant results
  go_output.SIG = go_output[go_output$over_represented_pvalue <= GOSEQ_SIGNIFICANCE_THRESHOLD | go_output$under_represented_pvalue <= GOSEQ_SIGNIFICANCE_THRESHOLD,]
  
  # Produce output table
  write.table(go_output.SIG, file=output.file.name, sep="\t", row.names = FALSE, quote=FALSE)
}
