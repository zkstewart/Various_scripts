# create_nmds.R
# Script to create a NMDS visualisation of read counts
# as would be used in a DGE analysis.

# The first thing you need to change are the two ALL-CAPS variables
# in the first "specify data file locations" section.

# You might need to change the row.names = bits in "data loading" section.

# You will need to add column names from your metadata to the
# data.scores table in the "NMDS" section under subheading
# "Add metadata from coldata to be used for plotting".




# Install packages
install.packages("ggplot2")
install.packages("vegan")

# Load packages
library(ggplot2)
library(vegan)

# Set working directory variable
WORK_DIR = "F:/this/can/be/anywhere"
setwd(WORK_DIR)


###########################################################################
##                                                                       ##
##                    SPECIFY DATA FILE LOCATIONS                        ##
##                                                                       ##
###########################################################################


COUNTS_FILE = "F:/this/can/be/anywhere/counts.tsv"
COLDATA_FILE = "F:/this/can/be/anywhere/coldata.txt"
## Remember to use forward slashes


###########################################################################
##                                                                       ##
##                           DATA LOADING                                ##
##                                                                       ##
###########################################################################


# Load in counts file
## If you've got a table like the below...
## +-----------------+------+------+-----+
## |     gene_id     | alm1 | alm2 | ... |
## +-----------------+------+------+-----+
## | Prudul26A002130 |  166 |  128 | ... |
## | Prudul26A026237 |  260 |  118 | ... |
## +-----------------+------+------+-----+
## ... then you should set row.names to be "gene_id"
counts.table = read.table(file=COUNTS_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = "gene_id")


# Load in coldata/metadata file
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



###########################################################################
##                                                                       ##
##                              NMDS                                     ##
##                                                                       ##
###########################################################################


# Transform table to work with vegan package NMDS
counts.table.t <- as.data.frame(t(counts.table))


# Perform the NMDS creation
## k=2 means we'll have two 2 dimensions/"principle components"
## trymax=100 should work unless you get a warning about solution convergence not occurring; in which case increase it
nmds=metaMDS(counts.table.t, k=2, trymax=100)
data.scores <- as.data.frame(scores(nmds)$sites)
data.scores$sample <- rownames(data.scores)


# Add metadata from coldata to be used for plotting
## Any column names from your coldata.table that you want to plot should be applied here
data.scores$tissue <- factor(coldata.table$tissue)
data.scores$time <- factor(coldata.table$time)
data.scores$replicate <- factor(coldata.table$replicate)


# Plot the NMDS
NMDS_FILE_NAME = "NMDS.pdf"
pdf(file=NMDS_FILE_NAME, width = 7, height = 7) # width and height are measured in inches; experiment with them
ggplot() + 
  geom_point(data=data.scores, aes(x=NMDS1, y=NMDS2, colour=time, shape=tissue), size=3) + # remove/change colour and/or shape
  #scale_shape_manual(values=1:nlevels(data.scores$time)) + # turn this ON if you need more than 6 shapes
  theme_bw()
dev.off()

