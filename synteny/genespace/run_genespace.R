# run_genespace.R
# Script to run the GENESPACE pipeline on files pre-formatted
# in the specified working directory

library(GENESPACE)

# Following guide at: https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/genespaceGuide.html


###########################################################################
##                                                                       ##
##                         ANALYSIS SETUP                                ##
##                                                                       ##
###########################################################################


# Specify working directory & MCScanX location
WORK_DIR = "/home/stewarz2/plant_group/glauca/genespace"
path2mcscanx = "/home/stewarz2/various_programs/MCScanX"


# Specify how many cores to run multi-threaded steps with
NCPUS = 12



###########################################################################
##                                                                       ##
##                          RUN GENESPACE                                ##
##                                                                       ##
###########################################################################


# Initialise GENESPACE parameters
gpar <- init_genespace(
  wd = WORK_DIR,
  nCores = NCPUS,
  diamondUltraSens = TRUE,
  path2mcscanx = path2mcscanx)


# Run GENESPACE pipeline
out <- run_genespace(gpar, overwrite = T)


# Save results
save.image(file="genespace_workspace.RData")
