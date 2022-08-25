# DESeq2_mangoCultivars_v1.R
# Script to run DESeq2 analysis for plant RNA-seq study
# of mango cultivars.

# Install packages
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("goseq")
# 
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("pheatmap")
# install.packages("RColorBrewer")
# install.packages("vegan")
# install.packages("comprehenr")
# install.packages("UpSetR")

# Load packages
library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(goseq)
library(RColorBrewer)
library(vegan)
library(grid)
library(gridExtra)
library(comprehenr)
library(UpSetR)


# Setup directory locations
# >> DGE-related variables
COUNTS_FILE = "F:/plant_group/mango_cultivars/counts/star/mango_cultivars_star_counts.tsv"
COLDATA_FILE = "F:/plant_group/mango_cultivars/metadata/coldata.txt"


# >> goseq-related variables
GOANNOT_FILE = "F:/plant_group/mango_cultivars/annotations/mango.fix.GO.mapping"
KEGGANNOT_FILE = "F:/plant_group/mango_cultivars/annotations/manindi_flc.aa_keggmap.tsv"
LENGTHS_FILE = "F:/plant_group/mango_cultivars/annotations/manindi_flc.gene.fix.lengths"


# >> Working directory variable
WORK_DIR = "F:/plant_group/mango_cultivars/DGE"
setwd(WORK_DIR)


# Specify significance threshold
SIGNIFICANCE_CUTOFF = 0.001
GOSEQ_SIGNIFICANCE_THRESHOLD = 0.05


# Specify clustering parameters
MIN_NUMBER_GENES = 5


# Load in data files for DGE analysis
# >> Counts file
counts.table = read.table(file=COUNTS_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = "gene_id")


# >> Coldata file
coldata.table = read.table(file=COLDATA_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = "SampleName")


# Specify as factor any relevant numeric variables that will be part of the experimental design
coldata.table$month <- factor(coldata.table$month)
coldata.table$replicate <- factor(coldata.table$replicate)
coldata.table$tc <- factor(coldata.table$tc)



###########################################################################
##                                                                       ##
##                     FUNCTION DECLARATIONS                             ##
##                                                                       ##
###########################################################################

subset_counts_and_coldata <- function(counts.table, coldata.table, subsetVariable, subsetValue)
{
  coldata.table.subset = coldata.table[coldata.table[[subsetVariable]] == subsetValue,]
  rowsToGrab = to_vec(
    for(rn in rownames(coldata.table.subset)) if(grepl("^[[:digit:]]+", rn)) paste0("X", rn) else(rn)
  )
  counts.table.subset = subset(counts.table, select=rowsToGrab)
  
  return(list("coldata" = coldata.table.subset, "counts" = counts.table.subset))
}

run_timecourse_dge <- function(counts.table, coldata.table, outputDirectory, outputFileName, designFormula, subsetVariable=NULL, subsetValue=NULL){
  ## designFormula should be formula(paste("~",x)) where x == the time variable's name
  ## subsetVariable, if specified, should be the column name for the variable to subset on
  ## subsetValue, if specified, should be the column value to subset on
  
  if (! is.null(subsetVariable))
  {
    subsetReturnList = subset_counts_and_coldata(counts.table, coldata.table, subsetVariable, subsetValue)
    coldata.table.tc = subsetReturnList$coldata
    counts.table.tc = subsetReturnList$counts
  } else
  {
    coldata.table.tc = coldata.table
    counts.table.tc = counts.table
  }
  
  # Create a DESeq2 object for timecourse analysis
  dds.tc <- DESeqDataSetFromMatrix(countData = counts.table.tc, colData = coldata.table.tc, design = designFormula)
  dds.tc <- estimateSizeFactors(dds.tc)
  
  # Drop low count genes
  nc.tc <- counts(dds.tc, normalized=TRUE)
  filter.tc <- rowSums(nc.tc >= 10) >= 2
  dds.tc <- dds.tc[filter.tc,]
  
  # Run DGE with timecourse LRT
  dds.tc <- DESeq(dds.tc, test = "LRT", reduced = ~1)
  res.tc <- results(dds.tc, alpha=SIGNIFICANCE_CUTOFF)
  res.tc <- res.tc[order(res.tc$padj),]
  res.tc <- na.omit(res.tc)
  
  # Get DE genes
  res.DE.tc <- res.tc[res.tc$padj < SIGNIFICANCE_CUTOFF,]
  
  # Output DEG table
  DE_FILE=paste0(outputDirectory, "/", outputFileName)
  write.table(as.data.frame(res.DE.tc), file=DE_FILE, sep="\t", quote=FALSE, col.names=NA)
  
  return(list("dds" = dds.tc, "res" = res.DE.tc))
}

run_pairwise_dge <- function(counts.table, coldata.table, designFormula, subsetVariable=NULL, subsetValue=NULL){
  ## designFormula should be formula(paste("~",x)) where x == the time variable's name
  ## subsetVariable, if specified, should be the column name for the variable to subset on
  ## subsetValue, if specified, should be the column value to subset on
  
  if (! is.null(subsetVariable))
  {
    subsetReturnList = subset_counts_and_coldata(counts.table, coldata.table, subsetVariable, subsetValue)
    coldata.table.pw = subsetReturnList$coldata
    counts.table.pw = subsetReturnList$counts
  } else
  {
    coldata.table.pw = coldata.table
    counts.table.pw = counts.table
  }
  
  # Create a DESeq2 object for pairwise analysis
  dds.pw <- DESeqDataSetFromMatrix(countData = counts.table.pw, colData = coldata.table.pw, design = designFormula)
  dds.pw <- estimateSizeFactors(dds.pw)
  
  # Drop low count genes
  nc.pw <- counts(dds.pw, normalized=TRUE)
  filter.pw <- rowSums(nc.pw >= 10) >= 2
  dds.pw <- dds.pw[filter.pw,]
  
  # Run DESeq2
  dds.pw <- DESeq(dds.pw)
  
  return(dds.pw)
}

write_pairwise_results_from_dds <- function(dds, combinations.table, contrastVariable, outputDirectory, fileSuffix, SIGNIFICANCE_CUTOFF) {
  # Extract pairwise comparisons
  for (i in 1:nrow(combinations.table))
  {
    # Extract relevant details from comparison
    group1 = combinations.table[i,][[1]]
    group2 = combinations.table[i,][[2]]
    
    # Get DE genes
    res.DE.pw <- get_pairwise_degs_from_dds(dds, contrastVariable, group1, group2, SIGNIFICANCE_CUTOFF)
    
    # Output DE genes table
    DE_FILE = paste0(outputDirectory, "/", group1, "_vs_", group2, fileSuffix)
    write.table(as.data.frame(res.DE.pw), file=DE_FILE, sep="\t", quote=FALSE, col.names=NA)
  }
}

get_pairwise_degs_from_dds <- function(dds, contrastVariable, group1, group2, SIGNIFICANCE_CUTOFF) {
  # Obtain pairwise comparison of interest
  res.pw <- results(dds, alpha=SIGNIFICANCE_CUTOFF, contrast=c(contrastVariable, group1, group2))
  
  # Clean up
  res.pw <- res.pw[order(res.pw$padj),]
  res.pw <- na.omit(res.pw)
  
  # Get DE genes
  res.DE.pw <- res.pw[res.pw$padj < SIGNIFICANCE_CUTOFF,]
  
  return(res.DE.pw)
}

get_pairwise_combinations <- function(coldata.table, variableName){
  
  pw.combinations = data.frame(var1=numeric(0), var2=numeric(0))
  
  firstVars = c()
  for (var1 in unique(coldata.table[[variableName]]))
  {
    firstVars = c(firstVars, var1)
    for (var2 in unique(coldata.table[[variableName]]))
    {
      # Skip self-comparison, and redundant comparison
      if (var2 %in% firstVars)
      {
        next
      }
      # Store result
      else
      {
        pw.combinations[nrow(pw.combinations)+1, ] = c(var1, var2)
      }
    }
  }
  
  return(pw.combinations)
}


run_timecourse_minmax <- function(de_ids, dds_object, coldata.table, timeVariable, HIGHCUTOFF=0.5, LOWCUTOFF=0.5){
  
  # Get a normalised counts table, subsetted to be specific to DEGs
  counts.table.norm = as.data.frame(counts(dds_object, normalized=TRUE))
  counts.table.norm = counts.table.norm[rownames(counts.table.norm) %in% de_ids,]
  
  # Get normalised (min-max) expression values
  counts.table.minmaxNorm = data.frame(gene=numeric(0), expression=numeric(0), time=numeric(0), pattern=character(0))
  for (i in 1:nrow(counts.table.norm))
  {
    row = counts.table.norm[i,]
    
    
    # Calculate median values per time
    times = c()
    for (time in unique(coldata.table[[timeVariable]]))
    {
      timeMedian = median(as.numeric(row[as.character(coldata.table[[timeVariable]]) == time]))
      times = c(times, timeMedian)
    }
    
    # Normalise to 0->1 range
    maxValue = max(times)
    minValue = min(times)
    
    normTimes = c()
    for (time in times)
    {
      normTime = (time - minValue) / (maxValue - minValue)
      normTimes = c(normTimes, normTime)
    }
    
    # Order the values
    orderedNormTimes = c()
    for (x in 1:length(unique(coldata.table[[timeVariable]])))
    {
      timeIndex = which(unique(coldata.table[[timeVariable]]) == x)
      normTime = normTimes[[timeIndex]]
      orderedNormTimes = c(orderedNormTimes, normTime)
    }
    
    # Get the expression pattern
    pattern = c()
    for (time in orderedNormTimes)
    {
      pattern = c(
        pattern,
        if (time >= HIGHCUTOFF) "U" else if (time <= LOWCUTOFF) "D" else "N"
      )
    }
    pattern = paste(pattern, collapse="-")
    
    # Store results in our table
    for (x in 1:length(unique(coldata.table[[timeVariable]])))
    {
      normTime = orderedNormTimes[[x]]
      counts.table.minmaxNorm[nrow(counts.table.minmaxNorm)+1, ] = c(rownames(row), normTime, x, pattern)
    }
  }
  return(counts.table.minmaxNorm)
}

curate_patterns_into_clusters <- function(counts.table.minmaxNorm, minimumMembers=15){
  clusterNum = 1
  counts.table.minmaxNorm$cluster = NA
  for (pattern in sort(unique(counts.table.minmaxNorm$pattern)))
  {
    # Count how many members share this pattern
    pattern.df = counts.table.minmaxNorm[counts.table.minmaxNorm$pattern == pattern,]
    patternMembers = length(unique(pattern.df$gene))
    
    # If it meets cutoff, set it as a cluster
    if (patternMembers >= minimumMembers)
    {
      counts.table.minmaxNorm[counts.table.minmaxNorm$pattern == pattern,]$cluster = clusterNum
      clusterNum = clusterNum + 1
    }
  }
  return(counts.table.minmaxNorm)
}

plot_clusters <- function(counts.table.minmaxNorm, outputDirectory, outputPrefix){
  for (cluster in unique(counts.table.minmaxNorm$cluster))
  {
    if (is.na(cluster)) {
      next
    }
    
    # Subset a dataframe for plotting
    cluster.df = na.omit(counts.table.minmaxNorm[counts.table.minmaxNorm$cluster == cluster,])
    cluster.df$expression = as.numeric(cluster.df$expression)
    cluster.df$time = as.factor(cluster.df$time)
    pattern = cluster.df$pattern[[1]] # all patterns will be the same, just get the first one
    
    # Create the plot
    p = ggplot(cluster.df, aes(x = time, y = expression)) + # group = gene if doing geom_line()
      geom_boxplot(color="firebrick3", size=1) +
      stat_summary(aes(group=1), fun=median, geom="line", colour="deepskyblue2", size=2) + 
      ggtitle(paste0("Pattern=", pattern, ";Genes=", length(unique(cluster.df$gene)))) +
      ylab("min-max normalised expression") +
      xlab("time (week num.)") +
      #scale_x_discrete(expand = c(0, 0.1)) + # removes padding at horizontal borders of plot
      theme_bw()
    
    # Save it
    CLUST_PLOT_FILE_NAME = paste0(outputDirectory, "/", outputPrefix, cluster, ".png")
    ggsave(file = CLUST_PLOT_FILE_NAME, p, width=7, height=5)
  }
}


run_goseq <- function(de_ids, geneList, length.vector, 
                      splitgoannot, splitkeggannot, GOSEQ_SIGNIFICANCE_THRESHOLD,
                      OUTPUT_DIR, GOSEQ_FILENAME, KEGGSEQ_FILENAME){
  
  # Create a vector identifying which genes are identified as DE in this analysis
  DE.vector <- geneList %in% de_ids
  DE.vector = DE.vector * 1 # converts TRUE to 1 and FALSE to 0
  DE.vector = as.integer(DE.vector)
  names(DE.vector) = geneList
  
  
  # Run PWF for length bias
  pwf <- nullp(DEgenes=DE.vector, bias.data=length.vector)
  
  
  # Run goseq for GOs and KEGGs
  go_output <- goseq(pwf, gene2cat = splitgoannot)
  kegg_output <- goseq(pwf, gene2cat = splitkeggannot)
  
  
  # Limit output to significant results
  go_output.SIG = go_output[go_output$over_represented_pvalue <= GOSEQ_SIGNIFICANCE_THRESHOLD | go_output$under_represented_pvalue <= GOSEQ_SIGNIFICANCE_THRESHOLD,]
  kegg_output.SIG = kegg_output[kegg_output$over_represented_pvalue <= GOSEQ_SIGNIFICANCE_THRESHOLD | kegg_output$under_represented_pvalue <= GOSEQ_SIGNIFICANCE_THRESHOLD,]
  
  
  # Produce output tables
  GOSEQ_FILE = paste0(OUTPUT_DIR, "/", GOSEQ_FILENAME)
  write.table(go_output.SIG, file=GOSEQ_FILE, sep="\t", row.names = FALSE, quote=FALSE)
  
  KEGGSEG_FILE = paste0(OUTPUT_DIR, "/", KEGGSEQ_FILENAME)
  write.table(kegg_output.SIG, file=KEGGSEG_FILE, sep="\t", row.names = FALSE, quote=FALSE)
}

run_goseq_on_patterns <- function(counts.table.minmaxNorm, geneList, length.vector, 
                                  splitgoannot, splitkeggannot, GOSEQ_SIGNIFICANCE_THRESHOLD,
                                  OUTPUT_DIR, FILE_PREFIX){
  # >> Clusters by number
  for (cluster in unique(counts.table.minmaxNorm$cluster))
  {
    if (is.na(cluster)) {
      next
    }
    
    # Get gene IDs for this cluster
    cluster.DE.ids = unique(na.omit(counts.table.minmaxNorm[counts.table.minmaxNorm$cluster == cluster,])$gene)
    
    # Run GOseq with these gene IDs
    run_goseq(de_ids=cluster.DE.ids, geneList, length.vector, 
              splitgoannot, splitkeggannot,
              GOSEQ_SIGNIFICANCE_THRESHOLD,
              OUTPUT_DIR=OUTPUT_DIR,
              GOSEQ_FILENAME=paste0(FILE_PREFIX, "_cluster_", cluster, "_GOseq.tsv"),
              KEGGSEQ_FILENAME=paste0(FILE_PREFIX, "_cluster_", cluster, "_KEGGseq.tsv"))
  }
  
  
  # >> Clusters by time and expression direction
  for (timeNum in 1:length(unique(counts.table.minmaxNorm$time)))
  {
    for (direction in c("U", "D"))
    {
      # Get all genes that fit this time/direction
      matching.df = counts.table.minmaxNorm[counts.table.minmaxNorm$time == timeNum & sapply(strsplit(as.character(counts.table.minmaxNorm$pattern), "\\-"), function(x) x[[timeNum]]) == direction, ]
      
      # Get gene IDs for this group
      matching.DE.ids = matching.df$gene
      
      # Run GOseq with these gene IDs
      run_goseq(de_ids=matching.DE.ids, geneList, length.vector, 
                splitgoannot, splitkeggannot,
                GOSEQ_SIGNIFICANCE_THRESHOLD,
                OUTPUT_DIR=OUTPUT_DIR,
                GOSEQ_FILENAME=paste0(FILE_PREFIX, "_time_", timeNum, "_direction_", direction, "_GOseq.tsv"),
                KEGGSEQ_FILENAME=paste0(FILE_PREFIX, "_time_", timeNum, "_direction_", direction, "_KEGGseq.tsv"))
    }
  }
}

###########################################################################
##                                                                       ##
##                        DATA EXPLORATION                               ##
##                                                                       ##
###########################################################################



# Create a general purpose DESeq2 object for exploration
dds <- DESeqDataSetFromMatrix(countData = counts.table, colData = coldata.table, design = ~ 1)
dds <- estimateSizeFactors(dds)


# Drop low count genes
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 2
dds <- dds[filter,]


#### Distances heatmap

# Calculate distances
vsd <- vst(dds, blind = FALSE) 
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )


# Create heatmap
rownames(sampleDistMatrix) <- paste( vsd$cultivar, vsd$tc, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors) # Need to save this manually in R studio with export



###########################################################################
##                                                                       ##
##                               QC                                      ##
##                                                                       ##
###########################################################################

## Data exploration does not reveal anything requiring removal. Hurrah!


###########################################################################
##                                                                       ##
##                              NMDS                                     ##
##                                                                       ##
###########################################################################

## Some code here credit to https://chrischizinski.github.io/rstats/vegan-ggplot2/


# Create output directory
NMDS_DIR = paste0(WORK_DIR, "/nmds")
dir.create(NMDS_DIR, showWarnings = FALSE)


# Transform table to work with vegan package NMDS
counts.table.t <- as.data.frame(t(counts.table))


# Perform the NMDS creation
nmds=metaMDS(counts.table.t,k=2,trymax=100)


# Extract relevant data for plotting
data.scores <- as.data.frame(scores(nmds)$sites)
data.scores$sample <- rownames(data.scores)
data.scores$cultivar <- coldata.table$cultivar
data.scores$tc <- coldata.table$tc
data.scores$group <- paste0(coldata.table$cultivar, "-", coldata.table$tc)


# Make variables factors if relevant
data.scores$sample <- factor(data.scores$sample)
data.scores$cultivar <- factor(data.scores$cultivar)
data.scores$tc <- factor(data.scores$tc)
data.scores$group <- factor(data.scores$group)

# Plot the NMDS
## Hull area plotting is not helpful; it will be skipped here

## Design a colour palette that's nice for colour blind folks (credit to http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/)
cbbPalette <- c("#000000", "#0072B2", "#009E73", "#F0E442", "#E69F00", "#D55E00", "#CC79A7")

## Plot it
NMDS_FILE_NAME=paste0(WORK_DIR, "/nmds/", "mac_NMDS_postQC.png")
png(filename=NMDS_FILE_NAME, width = 1280, height = 680)
ggplot() + 
  geom_point(data=data.scores, aes(x=NMDS1, y=NMDS2, colour=cultivar, shape=tc), size=3) + # shape=condition
  scale_shape_manual(values=1:nlevels(data.scores$tc)) +
  scale_fill_manual(values=cbbPalette) +
  scale_colour_manual(values=cbbPalette) +
  #geom_polygon(data=hull.data, aes(x=NMDS1, y=NMDS2, fill=cultivar, group=group), alpha=0.50) + # add the convex hulls
  coord_equal() +
  theme_bw()
dev.off()



###########################################################################
##                                                                       ##
##                  DIFFERENTIAL GENE EXPRESSION                         ##
##                                                                       ##
###########################################################################


# Load in the data again (just in case we've played around above)
save.image("cultivars_preDGE.RData")
#load("cultivars_preDGE.RData")



####### ANALYSIS 1 #######
####### Time-course for each cultivar #######

# Create output directory
TC_DE_DIR = paste0(WORK_DIR, "/timecourse_DE")
dir.create(TC_DE_DIR, showWarnings = FALSE)


# Get the design formula
designFormula.tc = formula(paste("~", "tc"))


# Run TC analyses
tc.dge.1.1 = run_timecourse_dge(counts.table, coldata.table, TC_DE_DIR, "12.052.037_tc_results.tsv", designFormula.tc, "cultivar", "12.052.037")
dds.tc.1.1 = tc.dge.1.1$dds
res.tc.1.1 = tc.dge.1.1$res

tc.dge.1.2 = run_timecourse_dge(counts.table, coldata.table, TC_DE_DIR, "Thai.wild_tc_results.tsv", designFormula.tc, "cultivar", "Thai.wild")
dds.tc.1.2 = tc.dge.1.2$dds
res.tc.1.2 = tc.dge.1.2$res

tc.dge.1.3 = run_timecourse_dge(counts.table, coldata.table, TC_DE_DIR, "Ampalam_tc_results.tsv", designFormula.tc, "cultivar", "Ampalam")
dds.tc.1.3 = tc.dge.1.3$dds
res.tc.1.3 = tc.dge.1.3$res

tc.dge.1.4 = run_timecourse_dge(counts.table, coldata.table, TC_DE_DIR, "Keitt_tc_results.tsv", designFormula.tc, "cultivar", "Keitt")
dds.tc.1.4 = tc.dge.1.4$dds
res.tc.1.4 = tc.dge.1.4$res

tc.dge.1.5 = run_timecourse_dge(counts.table, coldata.table, TC_DE_DIR, "1243_tc_results.tsv", designFormula.tc, "cultivar", "1243")
dds.tc.1.5 = tc.dge.1.5$dds
res.tc.1.5 = tc.dge.1.5$res

tc.dge.1.6 = run_timecourse_dge(counts.table, coldata.table, TC_DE_DIR, "Laurina.Lombok_tc_results.tsv", designFormula.tc, "cultivar", "Laurina.Lombok")
dds.tc.1.6 = tc.dge.1.6$dds
res.tc.1.6 = tc.dge.1.6$res


####### ANALYSIS 2 #######
####### Pairwise of each cultivar at each time point #######


# Create output directory
PW_DE_DIR = paste0(WORK_DIR, "/pairwise_DE")
dir.create(PW_DE_DIR, showWarnings = FALSE)


# Get the design formula
designFormula.pw = formula(paste("~", "cultivar"))


# Run pairwise analyses at each time point
dds.pw.1 = run_pairwise_dge(counts.table, coldata.table, designFormula.pw, "tc", "1")

dds.pw.2 = run_pairwise_dge(counts.table, coldata.table, designFormula.pw, "tc", "2")

dds.pw.3 = run_pairwise_dge(counts.table, coldata.table, designFormula.pw, "tc", "3")


# Get a table indicating all non-redundant pairwise combinations
pw.all.combinations = get_pairwise_combinations(coldata.table, "cultivar")


# Extract pairwise results to their own files
write_pairwise_results_from_dds(dds.pw.1, pw.all.combinations, "cultivar", PW_DE_DIR, "_time1_results.tsv", SIGNIFICANCE_CUTOFF)

write_pairwise_results_from_dds(dds.pw.2, pw.all.combinations, "cultivar", PW_DE_DIR, "_time2_results.tsv", SIGNIFICANCE_CUTOFF)

write_pairwise_results_from_dds(dds.pw.3, pw.all.combinations, "cultivar", PW_DE_DIR, "_time3_results.tsv", SIGNIFICANCE_CUTOFF)


###########################################################################
##                                                                       ##
##                          VENN DIAGRAM                                 ##
##                                                                       ##
###########################################################################

## For this project, we want to employ a Venn diagram-like process of
## getting genes that are commonly found as up or down-regulated across
## a timecourse in each cultivar.

## In this section, we'll do the set overlaps to determine what the Venn
## groups will be, and also make the plot.

## Note: when log2FoldChange is +ve, it's relative to the first sample
## specified in the contrast.


# Load in the data again (just in case we've played around above)
save.image("cultivars_preVenn.RData")
#load("cultivars_preVenn.RData")


# Create output directory
VENN_DE_DIR = paste0(WORK_DIR, "/venn_timecourse_DE")
dir.create(VENN_DE_DIR, showWarnings = FALSE)


# Get DE results for each cultivar
res.tc.1.1.DE = rownames(res.tc.1.1) # 12.052.037
res.tc.1.2.DE = rownames(res.tc.1.2) # Thai.wild
res.tc.1.3.DE = rownames(res.tc.1.3) # Ampalam
res.tc.1.4.DE = rownames(res.tc.1.4) # Keitt
res.tc.1.5.DE = rownames(res.tc.1.5) # 1243
res.tc.1.6.DE = rownames(res.tc.1.6) # Laurina.Lombok


# Format lists for use with Venn diagramming package
vennGenes = list(
  "12.052.037" = res.tc.1.1.DE,
  "Thai.wild" = res.tc.1.2.DE,
  "Ampalam" = res.tc.1.3.DE,
  "Keitt" = res.tc.1.4.DE,
  "1243" = res.tc.1.5.DE,
  "Laurina.Lombok" = res.tc.1.6.DE
)


# Extract relevant venn groupings for later GOseq
only_1.1 = setdiff(
  res.tc.1.1.DE, Reduce(union, list(res.tc.1.2.DE, res.tc.1.3.DE, res.tc.1.4.DE, res.tc.1.5.DE, res.tc.1.6.DE))
)

only_1.2 = setdiff(
  res.tc.1.2.DE, Reduce(union, list(res.tc.1.1.DE, res.tc.1.3.DE, res.tc.1.4.DE, res.tc.1.5.DE, res.tc.1.6.DE))
)

only_1.3 = setdiff(
  res.tc.1.3.DE, Reduce(union, list(res.tc.1.1.DE, res.tc.1.2.DE, res.tc.1.4.DE, res.tc.1.5.DE, res.tc.1.6.DE))
)

only_1.4 = setdiff(
  res.tc.1.4.DE, Reduce(union, list(res.tc.1.1.DE, res.tc.1.2.DE, res.tc.1.3.DE, res.tc.1.5.DE, res.tc.1.6.DE))
)

only_1.5 = setdiff(
  res.tc.1.5.DE, Reduce(union, list(res.tc.1.1.DE, res.tc.1.2.DE, res.tc.1.3.DE, res.tc.1.4.DE, res.tc.1.6.DE))
)

only_1.6 = setdiff(
  res.tc.1.6.DE, Reduce(union, list(res.tc.1.1.DE, res.tc.1.2.DE, res.tc.1.3.DE, res.tc.1.4.DE, res.tc.1.5.DE))
)

all_venn = Reduce(intersect, list(res.tc.1.1.DE, res.tc.1.2.DE, res.tc.1.3.DE, res.tc.1.4.DE, res.tc.1.5.DE, res.tc.1.6.DE))


# Plot as UpSet diagram
VENN_FILE_NAME=paste0(VENN_DE_DIR, "/timecourse_upset.pdf")
pdf(file=VENN_FILE_NAME, width = 12, height = 6)
upset(fromList(vennGenes), nsets=length(vennGenes), order.by = "freq")
dev.off()


# Write relevant venn groupings to individual files
only_1.1.table = res.tc.1.1[rownames(res.tc.1.1) %in% only_1.1,]
DE_1.1_FILE = paste0(VENN_DE_DIR, "/12.052.037_tc_only_results.tsv")
write.table(as.data.frame(only_1.1.table), file=DE_1.1_FILE, sep="\t", quote=FALSE, col.names=NA)

only_1.2.table = res.tc.1.2[rownames(res.tc.1.2) %in% only_1.2,]
DE_1.2_FILE = paste0(VENN_DE_DIR, "/Thai.wild_tc_only_results.tsv")
write.table(as.data.frame(only_1.2.table), file=DE_1.2_FILE, sep="\t", quote=FALSE, col.names=NA)

only_1.3.table = res.tc.1.3[rownames(res.tc.1.3) %in% only_1.3,]
DE_1.3_FILE = paste0(VENN_DE_DIR, "/Ampalam_tc_only_results.tsv")
write.table(as.data.frame(only_1.3.table), file=DE_1.3_FILE, sep="\t", quote=FALSE, col.names=NA)

only_1.4.table = res.tc.1.4[rownames(res.tc.1.4) %in% only_1.4,]
DE_1.4_FILE = paste0(VENN_DE_DIR, "/Keitt_tc_only_results.tsv")
write.table(as.data.frame(only_1.4.table), file=DE_1.4_FILE, sep="\t", quote=FALSE, col.names=NA)

only_1.5.table = res.tc.1.5[rownames(res.tc.1.5) %in% only_1.5,]
DE_1.5_FILE = paste0(VENN_DE_DIR, "/1243_tc_only_results.tsv")
write.table(as.data.frame(only_1.5.table), file=DE_1.5_FILE, sep="\t", quote=FALSE, col.names=NA)

only_1.6.table = res.tc.1.6[rownames(res.tc.1.6) %in% only_1.6,]
DE_1.6_FILE = paste0(VENN_DE_DIR, "/Laurina.Lombok_tc_only_results.tsv")
write.table(as.data.frame(only_1.6.table), file=DE_1.6_FILE, sep="\t", quote=FALSE, col.names=NA)

DE_all_venn_FILE = paste0(VENN_DE_DIR, "/all_tc_results.ids")
write.table(as.data.frame(all_venn), file=DE_all_venn_FILE, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)


###########################################################################
##                                                                       ##
##                             CLUSTERING                                ##
##                                                                       ##
###########################################################################


# Load in our environment if we've been modifying things elsewhere
save.image("cultivars_preClustering.RData")
#load("cultivars_preClustering.RData")


####### ANALYSIS 1 #######
####### Time-course for each cultivar #######


# Locate output directory
DE_CLUST_DIR = paste0(WORK_DIR, "/timecourse_DE/clustering")
dir.create(DE_CLUST_DIR, showWarnings = FALSE)


# Get normalised (min-max) expression values for each cultivar's timecourse DEGs
# >> 1.1
subsetResult1.1 = subset_counts_and_coldata(counts.table, coldata.table, "cultivar", "12.052.037")
counts.table.minmaxNorm.1.1 = run_timecourse_minmax(rownames(res.tc.1.1), dds.tc.1.1, subsetResult1.1$coldata, "tc")

# >> 1.2
subsetResult1.2 = subset_counts_and_coldata(counts.table, coldata.table, "cultivar", "Thai.wild")
counts.table.minmaxNorm.1.2 = run_timecourse_minmax(rownames(res.tc.1.2), dds.tc.1.2, subsetResult1.2$coldata, "tc")

# >> 1.3
subsetResult1.3 = subset_counts_and_coldata(counts.table, coldata.table, "cultivar", "Ampalam")
counts.table.minmaxNorm.1.3 = run_timecourse_minmax(rownames(res.tc.1.3), dds.tc.1.3, subsetResult1.3$coldata, "tc")

# >> 1.4
subsetResult1.4 = subset_counts_and_coldata(counts.table, coldata.table, "cultivar", "Keitt")
counts.table.minmaxNorm.1.4 = run_timecourse_minmax(rownames(res.tc.1.4), dds.tc.1.4, subsetResult1.4$coldata, "tc")

# >> 1.5
subsetResult1.5 = subset_counts_and_coldata(counts.table, coldata.table, "cultivar", "1243")
counts.table.minmaxNorm.1.5 = run_timecourse_minmax(rownames(res.tc.1.5), dds.tc.1.5, subsetResult1.5$coldata, "tc")

# >> 1.6
subsetResult1.6 = subset_counts_and_coldata(counts.table, coldata.table, "cultivar", "Laurina.Lombok")
counts.table.minmaxNorm.1.6 = run_timecourse_minmax(rownames(res.tc.1.6), dds.tc.1.6, subsetResult1.6$coldata, "tc")


# Assign cluster membership to genes
counts.table.minmaxNorm.1.1 = curate_patterns_into_clusters(counts.table.minmaxNorm.1.1)

counts.table.minmaxNorm.1.2 = curate_patterns_into_clusters(counts.table.minmaxNorm.1.2)

counts.table.minmaxNorm.1.3 = curate_patterns_into_clusters(counts.table.minmaxNorm.1.3)

counts.table.minmaxNorm.1.4 = curate_patterns_into_clusters(counts.table.minmaxNorm.1.4)

counts.table.minmaxNorm.1.5 = curate_patterns_into_clusters(counts.table.minmaxNorm.1.5)

counts.table.minmaxNorm.1.6 = curate_patterns_into_clusters(counts.table.minmaxNorm.1.6)


# Create plots for each cluster
plot_clusters(counts.table.minmaxNorm.1.1, DE_CLUST_DIR, "12.052.037_cluster_")

plot_clusters(counts.table.minmaxNorm.1.2, DE_CLUST_DIR, "Thai.wild_cluster_")

plot_clusters(counts.table.minmaxNorm.1.3, DE_CLUST_DIR, "Ampalam_cluster_")

plot_clusters(counts.table.minmaxNorm.1.4, DE_CLUST_DIR, "Keitt_cluster_")

plot_clusters(counts.table.minmaxNorm.1.5, DE_CLUST_DIR, "1243_cluster_")

plot_clusters(counts.table.minmaxNorm.1.6, DE_CLUST_DIR, "Laurina.Lombok_cluster_")


# Create output file
PATTERN_RESULTS_1.1_FILE_NAME = paste0(DE_CLUST_DIR, "/12.052.037_patterns_and_clusters.tsv")
write.table(counts.table.minmaxNorm.1.1, file=PATTERN_RESULTS_1.1_FILE_NAME, sep="\t", quote=FALSE, row.names=FALSE)

PATTERN_RESULTS_1.2_FILE_NAME = paste0(DE_CLUST_DIR, "/Thai.wild_patterns_and_clusters.tsv")
write.table(counts.table.minmaxNorm.1.2, file=PATTERN_RESULTS_1.2_FILE_NAME, sep="\t", quote=FALSE, row.names=FALSE)

PATTERN_RESULTS_1.3_FILE_NAME = paste0(DE_CLUST_DIR, "/Ampalam_patterns_and_clusters.tsv")
write.table(counts.table.minmaxNorm.1.3, file=PATTERN_RESULTS_1.3_FILE_NAME, sep="\t", quote=FALSE, row.names=FALSE)

PATTERN_RESULTS_1.4_FILE_NAME = paste0(DE_CLUST_DIR, "/Keitt_patterns_and_clusters.tsv")
write.table(counts.table.minmaxNorm.1.4, file=PATTERN_RESULTS_1.4_FILE_NAME, sep="\t", quote=FALSE, row.names=FALSE)

PATTERN_RESULTS_1.5_FILE_NAME = paste0(DE_CLUST_DIR, "/1243_patterns_and_clusters.tsv")
write.table(counts.table.minmaxNorm.1.5, file=PATTERN_RESULTS_1.5_FILE_NAME, sep="\t", quote=FALSE, row.names=FALSE)

PATTERN_RESULTS_1.6_FILE_NAME = paste0(DE_CLUST_DIR, "/Laurina.Lombok_patterns_and_clusters.tsv")
write.table(counts.table.minmaxNorm.1.6, file=PATTERN_RESULTS_1.6_FILE_NAME, sep="\t", quote=FALSE, row.names=FALSE)



###########################################################################
##                                                                       ##
##                              GOSEQ                                    ##
##                                                                       ##
###########################################################################


# Load in data in case we've been fiddling above (also, saves time for clustering)
save.image("cultivars_preGOSEQ.RData")
#load("cultivars_preGOSEQ.RData")


# Set up DEG vector
geneList = rownames(counts.table)


# Generate GO mapping structure
goannot=read.table(GOANNOT_FILE, as.is=T, header=F, sep='\t')
goannot[goannot == 0] <- NA
splitgoannot = strsplit(goannot[,2], split='; ')
name.list = as.vector(goannot[,1])
names(splitgoannot) = name.list
rm(goannot) # tidy up


# Generate KEGG mapping structure
keggannot=read.table(KEGGANNOT_FILE, as.is=T, header=F, sep='\t')
keggannot[keggannot == 0] <- NA
splitkeggannot = strsplit(keggannot[,2], split='; ')
name.list = as.vector(keggannot[,1])
names(splitkeggannot) = name.list
rm(keggannot) # tidy up


# Load in gene lengths as a vector
genelengths=read.table(LENGTHS_FILE, as.is=T, header=F)
prevector = genelengths[2]
rownames(prevector) = genelengths[,1] 
length.vector = as.integer(prevector[,1])
names(length.vector) = rownames(prevector)
rm(prevector) # tidy up
rm(genelengths) # tidy up


# Cull lengths that don't show up in the gene list
## This is relevant only for DEW which removes duplicate genes silently along the way
length.vector = length.vector[names(length.vector) %in% geneList]


####### ANALYSIS 1 #######
####### Time-course for each cultivar #######

run_goseq(de_ids=rownames(res.tc.1.1), geneList, length.vector, 
          splitgoannot, splitkeggannot,
          GOSEQ_SIGNIFICANCE_THRESHOLD,
          OUTPUT_DIR=paste0(WORK_DIR, "/timecourse_DE"),
          GOSEQ_FILENAME="12.052.037_tc_GOseq.tsv",
          KEGGSEQ_FILENAME="12.052.037_tc_KEGGseq.tsv")

run_goseq(de_ids=rownames(res.tc.1.2), geneList, length.vector, 
          splitgoannot, splitkeggannot,
          GOSEQ_SIGNIFICANCE_THRESHOLD,
          OUTPUT_DIR=paste0(WORK_DIR, "/timecourse_DE"),
          GOSEQ_FILENAME="Thai.wild_tc_GOseq.tsv",
          KEGGSEQ_FILENAME="Thai.wild_tc_KEGGseq.tsv")

run_goseq(de_ids=rownames(res.tc.1.3), geneList, length.vector, 
          splitgoannot, splitkeggannot,
          GOSEQ_SIGNIFICANCE_THRESHOLD,
          OUTPUT_DIR=paste0(WORK_DIR, "/timecourse_DE"),
          GOSEQ_FILENAME="Ampalam_tc_GOseq.tsv",
          KEGGSEQ_FILENAME="Ampalam_tc_KEGGseq.tsv")

run_goseq(de_ids=rownames(res.tc.1.4), geneList, length.vector, 
          splitgoannot, splitkeggannot,
          GOSEQ_SIGNIFICANCE_THRESHOLD,
          OUTPUT_DIR=paste0(WORK_DIR, "/timecourse_DE"),
          GOSEQ_FILENAME="Keitt_tc_GOseq.tsv",
          KEGGSEQ_FILENAME="Keitt_tc_KEGGseq.tsv")

run_goseq(de_ids=rownames(res.tc.1.5), geneList, length.vector, 
          splitgoannot, splitkeggannot,
          GOSEQ_SIGNIFICANCE_THRESHOLD,
          OUTPUT_DIR=paste0(WORK_DIR, "/timecourse_DE"),
          GOSEQ_FILENAME="1243_tc_GOseq.tsv",
          KEGGSEQ_FILENAME="1243_tc_KEGGseq.tsv")

run_goseq(de_ids=rownames(res.tc.1.6), geneList, length.vector, 
          splitgoannot, splitkeggannot,
          GOSEQ_SIGNIFICANCE_THRESHOLD,
          OUTPUT_DIR=paste0(WORK_DIR, "/timecourse_DE"),
          GOSEQ_FILENAME="Laurina.Lombok_tc_GOseq.tsv",
          KEGGSEQ_FILENAME="Laurina.Lombok_tc_KEGGseq.tsv")


####### ANALYSIS 2 #######
####### Pairwise of each cultivar at each time point #######


# Extract pairwise cultivar results
for (i in 1:nrow(pw.all.combinations))
{
  # Extract relevant details from comparison
  contrast = "cultivar"
  group1 = pw.all.combinations[i,][[1]]
  group2 = pw.all.combinations[i,][[2]]
  
  # Obtain pairwise DE genes for the three time points
  res.pw.1 <- get_pairwise_degs_from_dds(dds.pw.1, contrast, group1, group2, SIGNIFICANCE_CUTOFF)
  
  res.pw.2 <- get_pairwise_degs_from_dds(dds.pw.2, contrast, group1, group2, SIGNIFICANCE_CUTOFF)
  
  res.pw.3 <- get_pairwise_degs_from_dds(dds.pw.3, contrast, group1, group2, SIGNIFICANCE_CUTOFF)
  
  # Run the goseq for each time point
  run_goseq(de_ids=rownames(res.pw.1), geneList, length.vector, 
            splitgoannot, splitkeggannot,
            GOSEQ_SIGNIFICANCE_THRESHOLD,
            OUTPUT_DIR=paste0(WORK_DIR, "/pairwise_DE"),
            GOSEQ_FILENAME=paste0(group1, "_vs_", group2, "_time1_GOseq.tsv"),
            KEGGSEQ_FILENAME=paste0(group1, "_vs_", group2, "_time1_KEGGseq.tsv"))
  
  run_goseq(de_ids=rownames(res.pw.2), geneList, length.vector, 
            splitgoannot, splitkeggannot,
            GOSEQ_SIGNIFICANCE_THRESHOLD,
            OUTPUT_DIR=paste0(WORK_DIR, "/pairwise_DE"),
            GOSEQ_FILENAME=paste0(group1, "_vs_", group2, "_time2_GOseq.tsv"),
            KEGGSEQ_FILENAME=paste0(group1, "_vs_", group2, "_time2_KEGGseq.tsv"))
  
  run_goseq(de_ids=rownames(res.pw.3), geneList, length.vector, 
            splitgoannot, splitkeggannot,
            GOSEQ_SIGNIFICANCE_THRESHOLD,
            OUTPUT_DIR=paste0(WORK_DIR, "/pairwise_DE"),
            GOSEQ_FILENAME=paste0(group1, "_vs_", group2, "_time3_GOseq.tsv"),
            KEGGSEQ_FILENAME=paste0(group1, "_vs_", group2, "_time3_KEGGseq.tsv"))
}



####### ANALYSIS 3 #######
####### Venn groupings #######

run_goseq(de_ids=only_1.1, geneList, length.vector, 
          splitgoannot, splitkeggannot,
          GOSEQ_SIGNIFICANCE_THRESHOLD,
          OUTPUT_DIR=paste0(WORK_DIR, "/venn_timecourse_DE"),
          GOSEQ_FILENAME="12.052.037_tc_only_GOseq.tsv",
          KEGGSEQ_FILENAME="12.052.037_tc_only_KEGGseq.tsv")

run_goseq(de_ids=only_1.2, geneList, length.vector, 
          splitgoannot, splitkeggannot,
          GOSEQ_SIGNIFICANCE_THRESHOLD,
          OUTPUT_DIR=paste0(WORK_DIR, "/venn_timecourse_DE"),
          GOSEQ_FILENAME="Thai.wild_tc_only_GOseq.tsv",
          KEGGSEQ_FILENAME="Thai.wild_tc_only_KEGGseq.tsv")

run_goseq(de_ids=only_1.3, geneList, length.vector, 
          splitgoannot, splitkeggannot,
          GOSEQ_SIGNIFICANCE_THRESHOLD,
          OUTPUT_DIR=paste0(WORK_DIR, "/venn_timecourse_DE"),
          GOSEQ_FILENAME="Ampalam_tc_only_GOseq.tsv",
          KEGGSEQ_FILENAME="Ampalam_tc_only_KEGGseq.tsv")

run_goseq(de_ids=only_1.4, geneList, length.vector, 
          splitgoannot, splitkeggannot,
          GOSEQ_SIGNIFICANCE_THRESHOLD,
          OUTPUT_DIR=paste0(WORK_DIR, "/venn_timecourse_DE"),
          GOSEQ_FILENAME="Keitt_tc_only_GOseq.tsv",
          KEGGSEQ_FILENAME="Keitt_tc_only_KEGGseq.tsv")

run_goseq(de_ids=only_1.5, geneList, length.vector, 
          splitgoannot, splitkeggannot,
          GOSEQ_SIGNIFICANCE_THRESHOLD,
          OUTPUT_DIR=paste0(WORK_DIR, "/venn_timecourse_DE"),
          GOSEQ_FILENAME="1243_tc_only_GOseq.tsv",
          KEGGSEQ_FILENAME="1243_tc_only_KEGGseq.tsv")

run_goseq(de_ids=only_1.6, geneList, length.vector, 
          splitgoannot, splitkeggannot,
          GOSEQ_SIGNIFICANCE_THRESHOLD,
          OUTPUT_DIR=paste0(WORK_DIR, "/venn_timecourse_DE"),
          GOSEQ_FILENAME="Laurina.Lombok_tc_only_GOseq.tsv",
          KEGGSEQ_FILENAME="Laurina.Lombok_tc_only_KEGGseq.tsv")

run_goseq(de_ids=all_venn, geneList, length.vector, 
          splitgoannot, splitkeggannot,
          GOSEQ_SIGNIFICANCE_THRESHOLD,
          OUTPUT_DIR=paste0(WORK_DIR, "/venn_timecourse_DE"),
          GOSEQ_FILENAME="all_tc_only_GOseq.tsv",
          KEGGSEQ_FILENAME="all_tc_only_KEGGseq.tsv")


####### ANALYSIS 4 #######
####### Cluster groupings #######

# Locate output directory
GOSEQ_CLUST_DIR = paste0(WORK_DIR, "/timecourse_DE/clustering_goseq")
dir.create(GOSEQ_CLUST_DIR, showWarnings = FALSE)


# Run for each cluster and time direction on each cultivar's time course
run_goseq_on_patterns(counts.table.minmaxNorm.1.1,
                      geneList, length.vector, 
                      splitgoannot, splitkeggannot,
                      GOSEQ_SIGNIFICANCE_THRESHOLD,
                      OUTPUT_DIR=GOSEQ_CLUST_DIR,
                      FILE_PREFIX="12.052.037")

run_goseq_on_patterns(counts.table.minmaxNorm.1.2,
                      geneList, length.vector, 
                      splitgoannot, splitkeggannot,
                      GOSEQ_SIGNIFICANCE_THRESHOLD,
                      OUTPUT_DIR=GOSEQ_CLUST_DIR,
                      FILE_PREFIX="Thai.wild")

run_goseq_on_patterns(counts.table.minmaxNorm.1.3,
                      geneList, length.vector, 
                      splitgoannot, splitkeggannot,
                      GOSEQ_SIGNIFICANCE_THRESHOLD,
                      OUTPUT_DIR=GOSEQ_CLUST_DIR,
                      FILE_PREFIX="Ampalam")

run_goseq_on_patterns(counts.table.minmaxNorm.1.4,
                      geneList, length.vector, 
                      splitgoannot, splitkeggannot,
                      GOSEQ_SIGNIFICANCE_THRESHOLD,
                      OUTPUT_DIR=GOSEQ_CLUST_DIR,
                      FILE_PREFIX="Keitt")

run_goseq_on_patterns(counts.table.minmaxNorm.1.5,
                      geneList, length.vector, 
                      splitgoannot, splitkeggannot,
                      GOSEQ_SIGNIFICANCE_THRESHOLD,
                      OUTPUT_DIR=GOSEQ_CLUST_DIR,
                      FILE_PREFIX="1243")

run_goseq_on_patterns(counts.table.minmaxNorm.1.6,
                      geneList, length.vector, 
                      splitgoannot, splitkeggannot,
                      GOSEQ_SIGNIFICANCE_THRESHOLD,
                      OUTPUT_DIR=GOSEQ_CLUST_DIR,
                      FILE_PREFIX="Laurina.Lombok")



###########################################################################
##                                                                       ##
##                           FINAL NOTES                                 ##
##                                                                       ##
###########################################################################


# If you want further help, email me at zkstewart1@gmail.com.
