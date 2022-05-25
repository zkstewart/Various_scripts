# pca_plotter.R
# Script to visualise VCF SNP results

######
# APPROACH
# SNPRelate
# https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html#data-analysis
#####

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("gdsfmt")
# BiocManager::install("SNPRelate")

# Setup working directory
setwd("F:/flies/genome_based_2022/original/Parental_Selected/pca")

# Load library
library(SNPRelate)
library(ggplot2)
library(ggrepel)

# Specify input file locations
VCFFILE="F:/flies/genome_based_2022/original/Parental_Selected/freebayes/btrys06_parental.final.vcf"
SAMPLEIDSFILE="F:/flies/genome_based_2022/original/Parental_Selected/metadata/sampleids.txt"
SPECIESFILE="F:/flies/genome_based_2022/original/Parental_Selected/metadata/populations.txt"

# Convert VCF to workable format
GDSFILE="btrys06_parental.final.gds"
#snpgdsVCF2GDS(VCFFILE, GDSFILE, method="biallelic.only")
snpgdsSummary(GDSFILE)

# Load back in GDS file
genofile <- snpgdsOpen(GDSFILE)

# Load in metadata
sample.id <- scan(SAMPLEIDSFILE, what=character())
pop_code <- scan(SPECIESFILE, what=character())
head(cbind(sample.id, pop_code))

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

# Tabulate PCA data
## Can customise the pop=() value for different results
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

# Plot PCA (default)
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))

# Publication-quality plot (with customisation)
pdf("btrys06_parental_pca.raw.pdf", width = 6, height = 5)
ggplot(tab, aes(x = EV1, y = EV2, color = pop)) + ## Manual modification may be necessary
  geom_point(size = 3, alpha=0.5) + ggtitle("PCA - population") +
  labs(x=paste("EV1 (", toString(pca$eigenval[[1]]), "%)", sep=""), y=paste("EV2 (", toString(pca$eigenval[[2]]), "%)", sep="")) +
  theme_bw()
dev.off()

# Labels plot (to investigate specific points)
labelsPlot <- ggplot(tab, aes(x = EV1, y = EV2)) + geom_point() + geom_text_repel(aes(label = sample.id), size = 4)
labelsPlot
ggsave("btrys06_parental_pca.labels.pdf", width=15, height=15)
