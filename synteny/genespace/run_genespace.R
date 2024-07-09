# run_genespace.R
# Script to run the GENESPACE pipeline on files pre-formatted
# in the specified working directory

library(GENESPACE)

# Following guide at: https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/genespaceGuide.html


###########################################################################
##                                                                       ##
##                     FUNCTION DECLARATIONS                             ##
##                                                                       ##
###########################################################################

## See https://github.com/jtlovell/GENESPACE/blob/master/R/run_orthofinder.R
## Function re-defined herein so as to parallelise OrthoFinder's clustering step
## so it doesn't take so damn long to run


run_orthofinder <- function(gsParam, verbose = TRUE){
  ##############################################################################
  # 1. Get things set up
  # -- 1.1 combine genomeIDs and outgroups
  genomeIDs <- gsParam$genomeIDs
  if(!is.na(gsParam$outgroup)[1]){
    genomeIDs <- c(genomeIDs, gsParam$outgroup)
    genomeIDs <- genomeIDs[!duplicated(genomeIDs)]
  }
  
  
  # -- 1.2 get paths
  ofDir <- gsParam$paths$orthofinder
  tmpDir <- gsParam$paths$tmp
  pepDir <- gsParam$paths$peptide
  
  # -- 1.3 check to make sure that the peptide fastas are in there
  pepf <- file.path(pepDir, sprintf("%s.fa", genomeIDs))
  if(!all(file.exists(pepf)))
    stop(sprintf(
      "something is wrong with the peptide files. could not find: %s",
      paste(pepf[!file.exists(pepf)], collapse = "\n")))
  
  # -- 1.4 check that the orthofinder directory does not exist or is empty
  if(dir.exists(ofDir)){
    fs <- list.files(path = ofDir)
    if(length(fs) == 0){
      unlink(ofDir, recursive = T)
    }else{
      stop(sprintf(
        "orthofinder directory %s exists, remove this before proceeding\n",
        ofDir))
    }
  }
  
  # -- 1.5 rename the required parameters
  onewayBlast <- gsParam$params$onewayBlast
  diamondUltraSens <- gsParam$params$diamondUltraSens
  path2orthofinder <- gsParam$shellCalls$orthofinder
  nCores <- gsParam$params$nCores
  runOfInR <- !is.na(path2orthofinder)
  
  ############################################################################
  # 2. set up the directory structure and make sure things look good
  # -- copy peptides over to tmp directory
  # this allows for more genomeIDs in /peptide than just those in gsParam
  if(verbose)
    cat(strwrap(sprintf(
      "Copying files over to the temporary directory: %s",
      tmpDir), indent = 8, exdent = 16), sep = "\n")
  
  # -- if the tmpDir exists, remove and re-create it
  if(dir.exists(tmpDir))
    unlink(tmpDir, recursive = T)
  dir.create(tmpDir)
  
  # -- copy over the peptide files
  nu <- file.copy(pepf, tmpDir)
  
  ############################################################################
  # 3. Get the orthofinder command
  ofComm <- sprintf(
    "-f %s -t %s -a %s %s %s -X -o %s",
    tmpDir, nCores, nCores,
    ifelse(onewayBlast, "-1", ""),
    ifelse(diamondUltraSens, "-S diamond_ultra_sens", ""),
    ofDir)
  
  # -- strip out extra spaces that may exist
  ofComm <- gsub("  ", " ", gsub("  ", " ", ofComm))
  
  ############################################################################
  # 4. If orthofinder is available, run it from R
  if(runOfInR){
    if(verbose)
      cat(strwrap(sprintf(
        "Running the following command in the shell: `%s %s`.This can take a
        while. To check the progress, look in the `WorkingDirectory` in the
        output (-o) directory", path2orthofinder, ofComm),
        indent = 8, exdent = 16), sep = "\n")
    
    outp <- system2(
      path2orthofinder,
      ofComm,
      stdout = TRUE, stderr = TRUE)
    if(verbose)
      cat(paste(c("\t", outp), collapse = "\n\t"))
  }else{
    ############################################################################
    # 5. If not, print how to do it and stop.
    stop(cat(strwrap(
      "Could not find a valid path to the orthofinder program from R. To run
        orthofinder, ensure that the orthofinder program is in the $PATH, then
        call the following from the shell: \n", indent = 0, exdent = 8),
      sprintf("orthofinder %s", ofComm)),
      "Once OrthoFinder has been run, re-call run_genespace",sep = "\n")
  }
}



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
