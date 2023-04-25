# vcf_pca_plotter.R
# Script to read in a VCF and produce publication-quality
# PCA plots.


# Load packages
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("gdsfmt")
#BiocManager::install("SNPRelate")

#require("RColorBrewer")
library(SNPRelate)
library(ggplot2)
library(ggplot2)
library(ggrepel)



###########################################################################
##                                                                       ##
##                         ANALYSIS SETUP                                ##
##                                                                       ##
###########################################################################
## Refer to https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html#data-analysis
## for information on how to use SNPRelate


# Setup working directory
WORK_DIR="F:/daniel_mentorship/citrus/pca"
setwd(WORK_DIR)


# Specify file prefix for VCF / output files
PREFIX="citrus.final"


# Specify metadata file locations
METADATA_DIR = "F:/daniel_mentorship/citrus/metadata"
CITRUS_BULKS = paste0(METADATA_DIR, "/linkage_pops_full.txt") # 2 columns, SAMPLE_ID : BULK# format


# Specify VCF file location
VCFFILE=paste0("F:/daniel_mentorship/citrus/mpileup/", PREFIX, ".vcf")


# Load in metadata
citrus_bulks.table = read.table(file=CITRUS_BULKS, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(citrus_bulks.table) = c("code", "bulk")



###########################################################################
##                                                                       ##
##                          EXPLORE PCA                                  ##
##                                                                       ##
###########################################################################


# Convert VCF to workable format
GDSFILE=paste0(PREFIX, ".gds")
snpgdsVCF2GDS(VCFFILE, GDSFILE, method="biallelic.only")
snpgdsSummary(GDSFILE)


# Load back in GDS file
genofile <- snpgdsOpen(GDSFILE)


# LD-based SNP pruning
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, autosome.only=FALSE) # Try different LD thresholds for sensitivity analysis


# Get all selected snp id
snpset.id <- unlist(unname(snpset))
head(snpset.id)


# Run PCA
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only=FALSE)
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id")) # just get metadata now


# Assess PCA results
## variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))



###########################################################################
##                                                                       ##
##                      PUBLICATION PLOTTING                             ##
##                                                                       ##
###########################################################################


# Tabulate data
tab <- data.frame(sample.id = pca$sample.id,
                  bulk = citrus_bulks.table[match(pca$sample.id, citrus_bulks.table$code), ]$bulk,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)


# Write plot to PDF
PLOT_FILENAME=paste0(PREFIX, "_pca.raw.pdf")
pdf(PLOT_FILENAME, width = 6, height = 5)

ggplot(tab, aes(x = EV1, y = EV2, color = bulk)) + ## Manual modification may be necessary
  geom_point(size = 3, alpha=0.5) + ggtitle("PCA - bulks") +
  labs(x=paste("EV1 (", toString(pca$eigenval[[1]]), "%)", sep=""), y=paste("EV2 (", toString(pca$eigenval[[2]]), "%)", sep="")) +
  theme_bw()

dev.off()


# Write plot with outlier points labelled
OUTLIER_PLOT_FILENAME=paste0(PREFIX, "_pca.labelled.pdf")
pdf(OUTLIER_PLOT_FILENAME, width = 12, height = 10)

ggplot(tab, aes(x = EV1, y = EV2, color = bulk)) + ## Manual modification may be necessary
  geom_point(size = 3, alpha=0.5) + ggtitle("PCA - species") +
  labs(x=paste("EV1 (", toString(pca$eigenval[[1]]), "%)", sep=""), y=paste("EV2 (", toString(pca$eigenval[[2]]), "%)", sep="")) +
  theme_bw() +
  geom_text_repel(aes(label = sample.id), size = 4)

dev.off()
