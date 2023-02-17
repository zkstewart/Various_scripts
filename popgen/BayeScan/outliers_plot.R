# outliers_plot.R
# Script to produce figures indicating FST across contigs
# with outlier SNPs marked specifically

# Load libraries
library(ggplot2)


###########################################################################
##                                                                       ##
##                         ANALYSIS SETUP                                ##
##                                                                       ##
###########################################################################


# Setup file directories
WORK_DIR="F:/flies/chapa_2022/fst_per_site/individual"
BAYESCAN_DIR="F:/flies/chapa_2022/bayescan/individual"

setwd(WORK_DIR)


# Specify prefixes for input files
PREFIX="c2022"


# Get input file locations
## the VCF needs to be the original one updated to reflect the variants that
## made it through GESTE file creation; use get_bayescan_vcf.py from Various_scripts
metadataFile = "F:/flies/chapa_2022/metadata/population_map_species.txt"
vcfFile = "F:/flies/chapa_2022/freebayes/individual/c2022.final.bayescan.vcf"


# Auto-fill other file locations
fstFile = paste0(BAYESCAN_DIR, "/", PREFIX, "_fst.txt")
outlierFstFile = paste0(BAYESCAN_DIR, "/", PREFIX, ".outliers_fst.txt")



###########################################################################
##                                                                       ##
##                          DATA LOADING                                 ##
##                                                                       ##
###########################################################################


# Load in metadata
metadata = read.delim(metadataFile, header=FALSE)


# Load in BayeScan outputs
fstTable = read.delim(fstFile, header=TRUE, sep="")
outlierFstTable = read.delim(outlierFstFile, header=TRUE, sep="")

vcfTable = read.delim(vcfFile, header=FALSE, sep="\t", comment.char = '#')[, 1:2]
colnames(vcfTable) = c('CHROM', 'POS')


# Merge VCF details (chromosome, position) into FST table
fstTable$CHROM = vcfTable$CHROM
fstTable$POS = vcfTable$POS


# Merge outlier details into FST table
outliers = row.names(outlierFstTable)
fstTable$isOutlier = row.names(fstTable) %in% outliers


# Curate data: Retain only chromosomes with at least 10 SNPs
chromCount = table(fstTable$CHROM)
chromNames = names(chromCount)
keep <- vector(mode = "list")
CUTOFF=100
for (i in 1:length(chromCount)) {
  n = chromNames[[i]]
  c = chromCount[[i]]
  if (c > CUTOFF) {
    keep = c(keep, n)
  }
}
fstTable = subset(fstTable, CHROM %in% keep)


# Plot with default qplot stuff
qplot(POS, fst, facets = . ~ CHROM, col = factor(CHROM), data = fstTable) + theme_bw()


# Plot with customised ggplot2
p = ggplot(fstTable, aes(x = POS, y=fst, color=as.factor(isOutlier))) +
  geom_point() +
  facet_grid(rows=vars(CHROM), scales="free") +
  labs(x = "position", y = "FST") +
  theme_bw() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    legend.position = "none"
  )
p

# Save resulting plot
ggsave(paste0(PREFIX, ".fst_manhattan.png"), plot=last_plot(), height=6, width=10)

