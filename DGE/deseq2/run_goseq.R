# run_goseq.R
# Script to perform GOseq based on minimal required inputs
# The only section that requires modification is the 'ANALYSIS SETUP'
# portion. All other sections are intended to be left as-is.

library(DESeq2)
library(goseq)


###########################################################################
##                                                                       ##
##                     FUNCTION DECLARATIONS                             ##
##                                                                       ##
###########################################################################


run_goseq <- function(de_ids, geneList, length.vector, 
                      splitgoannot,
                      GOSEQ_SIGNIFICANCE_THRESHOLD){
  
  # Create a vector identifying which genes are identified as DE in this analysis
  DE.vector <- geneList %in% de_ids
  DE.vector = DE.vector * 1 # converts TRUE to 1 and FALSE to 0
  DE.vector = as.integer(DE.vector)
  names(DE.vector) = geneList
  
  # Run PWF for length bias
  pwf <- nullp(DEgenes=DE.vector, bias.data=length.vector)
  
  # Run goseq
  go_output <- goseq(pwf, gene2cat = splitgoannot)
  
  # Limit output to significant results
  go_output.SIG = go_output[go_output$over_represented_pvalue <= GOSEQ_SIGNIFICANCE_THRESHOLD | go_output$under_represented_pvalue <= GOSEQ_SIGNIFICANCE_THRESHOLD,]
  
  return(go_output.SIG)
}


###########################################################################
##                                                                       ##
##                         ANALYSIS SETUP                                ##
##                                                                       ##
###########################################################################


# Setup working directory variable
WORK_DIR = "F:/plant_group/james/GOseq"
setwd(WORK_DIR)


# Specify significance threshold
GOSEQ_SIGNIFICANCE_THRESHOLD = 0.05


# Load in GO annotations
## Table should be in a format like below:
## 1) no header row
## 2) two columns where genes (left column) lacking annotations
##    are denoted with '0', and annotation terms are separated with '; '
## +-------------+------------------------+
## |  cluster-0  | 0                      |
## |  cluster-1  | GO:0003954; GO:0110165 |
## |  cluster-2  | GO:0003954; ...        |
## +-------------+------------------------+
GOANNOT_FILE = "F:/plant_group/james/annotation/james_glauca_binge.GOs.tsv"


# Load in gene lengths
## Table should be in a format like below:
## 1) no header row
## 2) two columns where genes (left column) are matched up with
##    their CDS length
## +-------------+-------+
## |  cluster-0  |  386  |
## |  cluster-1  |  939  |
## |  cluster-2  |  420  |
## +-------------+-------+
LENGTHS_FILE = "F:/plant_group/james/annotation/BINge_clustered_result.onlyBinned.lengths_tsv"


# Load in gene IDs for testing
## Table should be in a format like below:
## 1) a header row where the first column is labelled 'ALL_IDS' and
##    lists every gene ID that was involved in the DGE analysis
## 2) one or more columns labelled with the experiment/test being
##    performed (any name is suitable), with the genes found to be
##    differentially expressed in that test listed below.
## +-------------+---------------+---------------+-----+
## |   ALL_IDS   |     test1     |     test2     | ... |
## +-------------+---------------+---------------+-----+
## |  cluster-0  |   cluster-0   |   cluster-7   | ... |
## |  cluster-1  |   cluster-42  |   cluster-33  | ... |
## |  cluster-2  |   cluster-76  |               | ... |
## |  cluster-3  |      ...      |               | ... |
## +-------------+---------------+---------------+-----+
IDS_FILE = "F:/plant_group/james/annotation/genes_to_test.tsv"



###########################################################################
##                                                                       ##
##                             GOSEQ SETUP                               ##
##                                                                       ##
###########################################################################


# Load in IDs for testing
ids.table = read.table(IDS_FILE, as.is=T, header=T, sep='\t')
geneList = ids.table$ALL_IDS


# Extract IDs for each individual test
tests.list = as.list(ids.table[, 2:ncol(ids.table),drop=FALSE])
tests.list = lapply(tests.list, function(x) x[x!=""])


# Generate GO mapping structure
goannot = read.table(GOANNOT_FILE, as.is=T, header=F, sep='\t')
goannot[goannot == 0] <- NA
splitgoannot = strsplit(goannot[,2], split='; ')
name.list = as.vector(goannot[,1])
names(splitgoannot) = name.list
rm(goannot) # tidy up


# Load in gene lengths as a vector
genelengths=read.table(LENGTHS_FILE, as.is=T, header=F)
prevector = genelengths[2]
rownames(prevector) = genelengths[,1] 
length.vector = as.integer(prevector[,1])
names(length.vector) = rownames(prevector)
rm(prevector) # tidy up
rm(genelengths) # tidy up


# Make sure geneList and all vectors are ordered equivalently
length.vector = length.vector[geneList]
splitgoannot = splitgoannot[geneList]
stopifnot(names(length.vector) == names(splitgoannot), names(length.vector) == geneList)



###########################################################################
##                                                                       ##
##                              RUN GOSEQ                                ##
##                                                                       ##
###########################################################################

for (testName in names(tests.list))
{
  goseq.result = run_goseq(tests.list[[testName]], geneList,
                           length.vector, splitgoannot,
                           GOSEQ_SIGNIFICANCE_THRESHOLD)
  write.table(goseq.result,
              file=paste0(testName, "_GOseq.tsv"),
              sep="\t", row.names = FALSE, quote=FALSE)
}



###########################################################################
##                                                                       ##
##                             DONE !!!                                  ##
##                                                                       ##
###########################################################################