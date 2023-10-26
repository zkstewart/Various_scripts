# pophelper_fs_display.R
# Script to visualise fastSTRUCTURE results

# Install and load library
remotes::install_github('royfrancis/pophelper')
library(pophelper)
library(grid)
library(gridExtra)
library(RColorBrewer)


###########################################################################
##                                                                       ##
##                         ANALYSIS SETUP                                ##
##                                                                       ##
###########################################################################


# Setup working directory
WORK_DIR = "F:/flies/chapa_2022/fastSTRUCTURE"
setwd(WORK_DIR)


# Specify directory containing fastStructure outputs (.meanQ file[s])
FS_OUT_DIR = "F:/flies/chapa_2022/fastSTRUCTURE/individual/simple"


# Specify prefix of output files, and K value determined by chooseK function
PREFIX = "c2022_simple"
K = 2


# Specify metadata file
## Example table format shown below
## sample_id values must correspond to sample IDs in your VCF / BED files.
## Also, make sure your first column does have the "sample_id" header.
## +-------------------+------+---------+-----+
## |     sample_id     |  env | species | ... |
## +-------------------+------+---------+-----+
## |    01BN_MT_502    |  MT  |   BN    | ... |
## |    01BN_MT_503    |  MT  |   BN    | ... |
## +-------------------+------+---------+-----+
METADATA_FILE = "F:/flies/chapa_2022/metadata/fs_metadata_table.txt"



###########################################################################
##                                                                       ##
##                          DATA LOADING                                 ##
##                                                                       ##
###########################################################################


# Parse in metadata
metadata.table = read.table(file=METADATA_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# Get the expected .meanQ file name and load it in
meanQFile = paste0(FS_OUT_DIR, "/", PREFIX, ".", K, ".meanQ")
qdata = readQ(meanQFile)
rownames(qdata[[1]]) = metadata.table$sample_id



###########################################################################
##                                                                       ##
##                     GENERATE STRUCTURE PLOTS                          ##
##                                                                       ##
###########################################################################


# Generate a colour palette with number of colours = K
colourPalette = brewer.pal(K, "Set1")[1:K] # [1:K] at the end handles the Warning message that occurs if K < 3
## try experimenting with other colour palettes by replacing "Set1" with others
## described at https://r-graph-gallery.com/38-rcolorbrewers-palettes.html


# Specify the title for your plot
plotTitle = "fastSTRUCTURE admixture proportions for BNEO and BTRY"


# Generate simple plot of just sample IDs
plotQ(qlist=qdata, clustercol=colourPalette,
      # Chart size
      width = 12, height = 3, # try to experiment with sizes
      
      # Barplot aesthetics
      barbordersize = 0.1, barbordercolour="black",
      
      # Title
      showtitle = TRUE,
      titlelab = plotTitle,
      
      # Individual labels
      showindlab=T, useindlab=T, indlabsize=2.5,
      indlabheight=1, indlabvjust=1,
      
      # Sorting
      sortind="Cluster1", # remove if you want to use the default order of samples
      
      # Output plot
      imgtype="pdf",
      exportpath=WORK_DIR)


# Generate plot with additional labeling
plotQ(qlist=qdata, clustercol=colourPalette,
      # Chart size
      width = 18, height = 3, # try to experiment with sizes
      
      # Barplot aesthetics
      barbordersize = 0.1, barbordercolour="black",
      
      # Title
      showtitle = TRUE,
      titlelab = plotTitle,
      
      # Individual labels
      showindlab=T, useindlab=T, indlabsize=2.5,
      indlabheight=1, indlabvjust=1,
      
      # Group labels
      grplab=metadata.table$environment, # change this to whatever column you want from your metadata!
      grplabsize=2, linesize=0.8, pointsize=3,
      grplabangle=90, grplabheight=2.5,
      
      # Sorting
      ordergrp=T, # should probably use this if plotting with extra metadata
      
      # Output plot
      imgtype="pdf",
      exportpath=WORK_DIR)


###########################################################################
##                                                                       ##
##                  OPTIONAL: PLOT K VALUE COMPARISON                    ##
##                                                                       ##
###########################################################################
## No need to use this code unless you're trying to decide between more than
## one K value. This plot should help to visualise the differences.


# Specify the other K value's .meanQ file
OTHER_K = 3

otherMeanQFile = paste0(FS_OUT_DIR, "/", PREFIX, ".", OTHER_K, ".meanQ")
otherQdata = readQ(otherMeanQFile)
rownames(otherQdata[[1]]) = metadata.table$sample_id


# Load multiple Q value files in together
fsFiles <- c(meanQFile, otherMeanQFile)
fsList <- readQ(files=fsFiles)
fsList <- lapply(fsList,"rownames<-",metadata.table$sample_id)


# Align the different K values
alignList <- alignK(fsList)


# Generate a colour palette with number of colours = the highest K
colourPalette = brewer.pal(max(K, OTHER_K), "Set1")[1:max(K, OTHER_K)]


# Plot both K values
plotQ(qlist=alignList, clustercol=colourPalette,
      # Chart size
      width = 12, height = 3, # try to experiment with sizes
      
      # Barplot aesthetics
      barbordersize = 0.1, barbordercolour="black",
      
      # Title
      showtitle = TRUE,
      titlelab = plotTitle,
      
      # Individual labels
      showindlab=T, useindlab=T, indlabsize=2.5,
      indlabheight=1, indlabvjust=1,
      
      # Output plot
      imgtype="pdf", imgoutput="join",
      exportpath=getwd())

