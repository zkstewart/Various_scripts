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
WORK_DIR="F:/flies/chapa_2022/pca/individual"
setwd(WORK_DIR)


# Specify file prefix for VCF / output files
PREFIX="c2022.final"


# Specify metadata file locations
SAMPLEIDSFILE="F:/flies/chapa_2022/metadata/fs_sampleids.txt" # 1 column, SAMPLE_CODE format
SPECIESFILE="F:/flies/chapa_2022/metadata/fs_species.txt" # 1 column, SP format
ENVIRONFILE="F:/flies/chapa_2022/metadata/fs_environ.txt" # 1 column, ENV format
SPENVFILE="F:/flies/chapa_2022/metadata/fs_spenv.txt" # 1 column, SP_ENV format


# Specify VCF file location
VCFFILE=paste0("F:/flies/chapa_2022/freebayes/individual/", PREFIX, ".vcf")


# Load in metadata
sample.id <- scan(SAMPLEIDSFILE, what=character())
sp_pop_code <- scan(SPECIESFILE, what=character())
env_pop_code <- scan(ENVIRONFILE, what=character())
spenv_pop_code <- scan(SPENVFILE, what=character())



###########################################################################
##                                                                       ##
##                          EXPLORE PCA                                  ##
##                                                                       ##
###########################################################################


# Convert VCF to workable format
GDSFILE="btrys06.final.gds"
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


# Plot PCA (species)
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(sp_pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomleft", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))


# Plot PCA (environ)
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(env_pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomleft", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))


# Plot PCA (sp_env)
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(spenv_pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomleft", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))


###########################################################################
##                                                                       ##
##                      PUBLICATION PLOTTING                             ##
##                                                                       ##
###########################################################################


# Tabulate data
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(sp_pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)


# Write plot to PDF
PLOT_FILENAME=paste0(PREFIX, "_pca.raw.pdf")
pdf(PLOT_FILENAME, width = 6, height = 5)

ggplot(tab, aes(x = EV1, y = EV2, color = pop)) + ## Manual modification may be necessary
  geom_point(size = 3, alpha=0.5) + ggtitle("PCA - species") +
  labs(x=paste("EV1 (", toString(pca$eigenval[[1]]), "%)", sep=""), y=paste("EV2 (", toString(pca$eigenval[[2]]), "%)", sep="")) +
  theme_bw()

dev.off()


# Write plot with outlier points labelled
OUTLIER_PLOT_FILENAME=paste0(PREFIX, "_pca.labelled.pdf")
pdf(OUTLIER_PLOT_FILENAME, width = 6, height = 5)

ggplot(tab, aes(x = EV1, y = EV2, color = pop)) + ## Manual modification may be necessary
  geom_point(size = 3, alpha=0.5) + ggtitle("PCA - species") +
  labs(x=paste("EV1 (", toString(pca$eigenval[[1]]), "%)", sep=""), y=paste("EV2 (", toString(pca$eigenval[[2]]), "%)", sep="")) +
  theme_bw() +
  geom_text_repel(aes(label = sample.id), size = 4)

dev.off()
