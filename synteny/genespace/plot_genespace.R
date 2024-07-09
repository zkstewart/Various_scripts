# plot_genespace.R
# Script to take the RData produced from a GENESPACE run
# and generate customised riparian plots

library(GENESPACE)
library(RColorBrewer)
library(ggplot2)

# Following guide at: https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/riparianGuide.html


###########################################################################
##                                                                       ##
##                         ANALYSIS SETUP                                ##
##                                                                       ##
###########################################################################

## Download all the folders from the GENESPACE run except for the tmp
## folder and the OrthoFinder subfolders with heaps of files
## and place them in WORK_DIR


# Setup working directory variable
WORK_DIR = "F:/plant_group/glauca/genespace/all"
setwd(WORK_DIR)


# Load workspace
WORKSPACE_FILE = "F:/plant_group/glauca/genespace/all/genespace_workspace.RData"
load(WORKSPACE_FILE)


# Modify paths to point to local file structure
OLD_DIR = "/home/stewarz2/plant_group/glauca/genespace" # all OLD_DIR should become WORK_DIR

# $paths
for (subname in names(out$paths))
{
  out$paths[[subname]] = gsub(OLD_DIR, WORK_DIR, out$paths[[subname]])
}

# $synteny$blast
for (columnname in c("queryBlast", "targetBlast", "queryOrthologs", "targetOrthologs", "allBlast", "synHits"))
{
  out$synteny$blast[[columnname]] = ifelse(is.na(out$synteny$blast[[columnname]]),
                                           NA,
                                           gsub(OLD_DIR, WORK_DIR, out$synteny$blast[[columnname]]))
}

# $synteny subnames
for (subname in c("combBed", "SpeciesIDs", "SequenceIDs", "ogs", "hogs", "speciesTree"))
{
  out$synteny[[subname]] = gsub(OLD_DIR, WORK_DIR, out$synteny[[subname]])
}


# Prevent Windows-related issues
out$params$nCores = 1


###########################################################################
##                                                                       ##
##                         RIPARIAN PLOTS                                ##
##                                                                       ##
###########################################################################


# Create directory for new output plots
PLOTS_DIR = paste0(WORK_DIR, "/custom_plots")
dir.create(PLOTS_DIR, showWarnings = FALSE)


# Get a colour palette with 9 values (for the 9 chromosomes)
myPalette = brewer.pal(12, name="Paired")[c(1:2,5:10, 12)] # cut out the green tones and the light yellow


# Generate plots
## Aussie only
ripd <- plot_riparian(
  gsParam = out,
  refGenome = "glauca",
  useRegions = FALSE,
  genomeIDs = c("australasica", "australis", "glauca"),
  inversionColor = "green",
  chrFill = "lightgrey",
  palette = colorRampPalette(myPalette),
  addThemes = ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white")),
  minChrLen2plot = 2000)
ggsave(filename=paste0(PLOTS_DIR, "/riparian_plot.pdf"),
       width=12, height=5)


## CAFE matched
invchr <- data.frame(
  genome = c(rep("hindsii", 4)), 
  chr = c(c(1, 3, 6, 7)))


ripd <- plot_riparian(
  gsParam = out,
  refGenome = "glauca",
  useRegions = FALSE,
  genomeIDs = c("sinensis", "limon", "hindsii", "australasica", "australis", "glauca"),
  invertTheseChrs = invchr,
  inversionColor = "green",
  chrFill = "lightgrey",
  palette = colorRampPalette(myPalette),
  addThemes = ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white")),
  minChrLen2plot = 2000)
ggsave(filename=paste0(PLOTS_DIR, "/riparian_cafespecies_plot.pdf"),
       width=12, height=8)

