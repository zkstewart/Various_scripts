# plot_faststructure.R
# Script to visualise fastSTRUCTURE results

# Install and load library
# remotes::install_github('royfrancis/pophelper')
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
WORK_DIR = "F:/lab_members/devindee/structure"
setwd(WORK_DIR)

# Specify directory containing fastStructure outputs (.meanQ file[s])
FS_OUT_DIR = "F:/lab_members/devindee/structure"

# Specify prefix of output files, and K value determined by chooseK function
PREFIX = "flgenomics_commercial_simple"
OPTIMAL_K = 2

# Specify .fam metadata file
## This file is likely to have been used during previous steps converting
## your VCF to PLINK2 BED format.
FAM_FILE = "F:/lab_members/devindee/metadata/flgenomics_commercial.fam"

# Load in .fam metadata
metadata.table = read.table(file=FAM_FILE, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(metadata.table) = c("FID", "IID", "PARENT1", "PARENT2", "SEX", "PHENO")


###########################################################################
##                                                                       ##
##                     GENERATE STRUCTURE PLOTS                          ##
##                                                                       ##
###########################################################################

# Generate a colour palette with number of colours = K
## try experimenting with other colour palettes by replacing "Set1" with others
## described at https://r-graph-gallery.com/38-rcolorbrewers-palettes.html
colourPalette = brewer.pal(OPTIMAL_K, "Set1")[1:OPTIMAL_K] # [1:OPTIMAL_K] at the end handles the Warning message that occurs if K < 3

# Get the optimal K .meanQ file name and load it in
meanQFile = file.path(FS_OUT_DIR, paste0(PREFIX, ".", OPTIMAL_K, ".meanQ"))
qdata = readQ(meanQFile)
rownames(qdata[[1]]) = metadata.table$IID

# Specify the title for your plot
plotTitle = "Admixture proportions for RL and commercial FL varieties"

# Generate simple plot of just sample IDs
plotQ(qlist=qdata, clustercol=colourPalette,
      # Chart size
      width = 12, height = 3, # try to experiment with sizes
      
      # Barplot aesthetics
      barbordersize = 0.1,
      barbordercolour="black",
      
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
      barbordersize = 0.1,
      barbordercolour="black",
      
      # Title
      showtitle = TRUE,
      titlelab = plotTitle,
      
      # Individual labels
      showindlab=T, useindlab=T, indlabsize=2.5,
      indlabheight=1, indlabvjust=1,
      
      # Group labels
      grplab=data.frame(metadata.table$FID), # change this to whatever column you want from your metadata!
      grplabsize=1.5, linesize=0.8, pointsize=3,
      grplabangle=90, grplabheight=2.5,
      grplabjust=0.25,
      
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

otherMeanQFile = file.path(FS_OUT_DIR, paste0(PREFIX, ".", OTHER_K, ".meanQ"))
otherQdata = readQ(otherMeanQFile)
rownames(otherQdata[[1]]) = metadata.table$IID

# Load multiple Q value files in together
fsFiles <- c(meanQFile, otherMeanQFile)
fsList <- readQ(files=fsFiles)
fsList <- lapply(fsList,"rownames<-",metadata.table$IID)

# Align the different K values
alignList <- alignK(fsList)

# Generate a colour palette with number of colours = the highest K
colourPalette = brewer.pal(max(OPTIMAL_K, OTHER_K), "Set1")[1:max(OPTIMAL_K, OTHER_K)]

# Plot both K values
plotQ(qlist=alignList, clustercol=colourPalette,
      # Chart size
      width = 12, height = 3, # try to experiment with sizes
      
      # Barplot aesthetics
      barbordersize = 0.1,
      barbordercolour="black",
      
      # Title
      showtitle = TRUE,
      titlelab = plotTitle,
      
      # Individual labels
      showindlab=T, useindlab=T, indlabsize=2.5,
      indlabheight=1, indlabvjust=1,
      
      # Output plot
      imgtype="pdf", imgoutput="join",
      exportpath=getwd())

# Plot both K values with additional labeling
plotQ(qlist=alignList, clustercol=colourPalette,
      # Chart size
      width = 18, height = 3, # try to experiment with sizes
      
      # Barplot aesthetics
      barbordersize = 0.1,
      barbordercolour="black",
      
      # Title
      showtitle = TRUE,
      titlelab = plotTitle,
      
      # Individual labels
      showindlab=T, useindlab=T, indlabsize=2.5,
      indlabheight=1, indlabvjust=1,
      
      # Group labels
      grplab=data.frame(metadata.table$FID), # change this to whatever column you want from your metadata!
      grplabsize=1.5, linesize=0.8, pointsize=3,
      grplabangle=90, grplabheight=2.5,
      grplabjust=0.25,
      
      # Sorting
      ordergrp=T, # should probably use this if plotting with extra metadata
      
      # Output plot
      imgtype="pdf", imgoutput="join",
      exportpath=WORK_DIR)
