library(ggplot2)
library(ggrepel)

###########################################################################
##                                                                       ##
##                         ANALYSIS SETUP                                ##
##                                                                       ##
###########################################################################

# Set working directory
WORK_DIR="F:/lab_members/devindee/pca/testrun2"
setwd(WORK_DIR)

# Specify the PCA result files
SSCORE_FILE = "flgenomics_commercial_projected.sscore"
EIGENVAL_FILE = "flgenomics_commercial_counts.eigenval"

# Specify prefix for output files
PREFIX = "flgenomics_commercial"

###########################################################################
##                                                                       ##
##                            PLOT PCA                                   ##
##                                                                       ##
###########################################################################

# Load result files
sscore.table = read.table(file=SSCORE_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
eigenval.list = read.table(file=EIGENVAL_FILE, header = FALSE, stringsAsFactors = FALSE)$V1

# Interpret eigenvalues as proportions of explained variance
total = sum(eigenval.list)
explained = round((eigenval.list / total) * 100, digits=2)

# Tabulate data
tab <- data.frame(sample.id = sscore.table$IID,
                  Family = sscore.table$X.FID,
                  PC1 = sscore.table$PC1_AVG,
                  PC2 = sscore.table$PC2_AVG,
                  stringsAsFactors = FALSE)

# Write plot to PDF
PLOT_FILENAME=paste0(PREFIX, ".rplot.pdf")
pdf(PLOT_FILENAME, width = 8, height = 7)

ggplot(tab, aes(x = PC1, y = PC2, color = Family)) +
  geom_point(size = 3, alpha=0.6) + ggtitle("PCA - Family") +
  labs(x=paste("PC1 (", toString(explained[[1]]), "%)", sep=""), y=paste("PC2 (", toString(explained[[2]]), "%)", sep="")) +
  theme_bw()

dev.off()
