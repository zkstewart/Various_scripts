# pophelper_fs_display.R
# Script to visualise fastSTRUCTURE results

# Install and load library
remotes::install_github('royfrancis/pophelper')
library(pophelper)
library(grid)
library(gridExtra)

# Setup working directory
setwd("F:/flies/chapa_2022/fastSTRUCTURE")

###### 
# K=2 ONLY PLOT
######

# Specify input file locations
meanQFile="btrys06_simple.2.meanQ"
labelsFile="F:/flies/chapa_2022/metadata/fs_sampleids.txt"
labelsFile_2="F:/flies/chapa_2022/metadata/fs_2labels.txt"

# Read in Q values
qdata = readQ(meanQFile)

# Read in individual labels & associate to qdata
pops = read.delim(labelsFile, header=FALSE)

## Alt presentation: only genotype number
newpops = c()
for (i in 1:length(rownames(pops))) {
  oldlabel = pops$V1[i]
  newlabel = unlist(strsplit(toString(oldlabel), split="_")[1])[3]
  if (newlabel == "E") {
    newlabel = unlist(strsplit(toString(oldlabel), split="_")[1])[4]
  }
  newpops = c(newpops, newlabel)
}
pops$V1 = newpops

rownames(qdata[[1]]) = pops$V1
rownames(qdata[["btrys06_simple.2"]]) = pops$V1

# Plot it
plotQ(qlist=qdata, clustercol=c("slateblue1", "firebrick3"),
      # Chart size
      width=12, height=3,
      # Barplot aesthetics
      barbordersize=0.1, barbordercolour="black",
      # Title
      showtitle=TRUE,
      titlelab="fastSTRUCTURE admixture proportions for BNEO and BTRY",
      # Individual labels
      showindlab=T, useindlab=T, indlabsize=2.5,
      indlabheight=1, indlabvjust=1,
      # Output type
      imgtype="pdf",
      # Necessary for export
      exportpath=getwd())

# Sorted plot
plotQ(qlist=qdata, clustercol=c("slateblue1", "firebrick3"),
      # Chart size
      width=12, height=3,
      # Barplot aesthetics
      barbordersize=0.1, barbordercolour="black",
      # Title
      showtitle=TRUE,
      titlelab="fastSTRUCTURE admixture proportions for BNEO and BTRY",
      # Individual labels
      showindlab=T, useindlab=T, indlabsize=2.5,
      indlabheight=1, indlabvjust=1,
      # Output type
      imgtype="pdf",
      # Sorting
      sortind="Cluster1",
      # Necessary for export
      exportpath=getwd())

# Read in group labels
labels2 = read.delim(labelsFile_2, header=FALSE, stringsAsFactors=F)
colnames(labels2) = c("environ", "species")
environLabels <- labels2[,1,drop=FALSE]

# Sorted, environment label plot
plotQ(qlist=qdata, clustercol=c("slateblue1", "firebrick3"),
      # Chart size
      width=18, height=3,
      # Barplot aesthetics
      barbordersize=0.1, barbordercolour="black",
      # Title
      showtitle=TRUE,
      titlelab="fastSTRUCTURE admixture proportions for BNEO and BTRY",
      # Individual labels
      showindlab=T, useindlab=T, indlabsize=2.5,
      indlabheight=1, indlabvjust=1,
      # Output type
      imgtype="pdf",
      # Group labels
      grplab=environLabels,grplabsize=2,linesize=0.8,pointsize=3,
      grplabangle=90, grplabheight=2.5,
      # Sorting
      ordergrp=T,
      # Necessary for export
      exportpath=getwd())

# Sorted, multigroup label plot
plotQ(qlist=qdata, clustercol=c("slateblue1", "firebrick3"),
      # Chart size
      width=18, height=3,
      # Barplot aesthetics
      barbordersize=0.1, barbordercolour="black",
      # Title
      showtitle=TRUE,
      titlelab="fastSTRUCTURE admixture proportions for BNEO and BTRY",
      # Individual labels
      showindlab=T, useindlab=T, indlabsize=2.5,
      indlabheight=1, indlabvjust=1,
      # Output type
      imgtype="pdf",
      # Group labels
      grplab=labels2,grplabsize=2,linesize=0.8,pointsize=3,
      grplabangle=90, grplabheight=2.5,
      # Sorting
      ordergrp=T,
      # Necessary for export
      exportpath=getwd())


###### 
# K=2 VS K=3 PLOT
######

# Specify input file locations
labelsFile="F:/flies/chapa_2022/metadata/fs_sampleids.txt"
k2File="btrys06_simple.2.meanQ"
k3File="btrys06_simple.3.meanQ"
fsFiles <- c(k2File, k3File)

# Read in Q values
fsList <- readQ(files=fsFiles)

# Read in individual labels & associate to qdata
pops = read.delim(labelsFile, header=FALSE)
rownames(fsList[[1]]) = pops$V1
if(length(unique(sapply(fsList,nrow)))==1) fsList <- lapply(fsList,"rownames<-",pops$V1)

# Align K=2 and K=3 runs
alignList <- alignK(fsList)

# Plot it
plotQ(qlist=alignList, clustercol=c("slateblue1", "firebrick3", "green"),
      # Chart size
      width=12, height=3,
      # Barplot aesthetics
      barbordersize=0.1, barbordercolour="black",
      # Title
      showtitle=TRUE,
      titlelab="fastSTRUCTURE admixture proportions for BNEO and BTRY",
      # Individual labels
      showindlab=T, useindlab=T, indlabsize=2.5,
      indlabheight=1, indlabvjust=1,
      # Output type
      imgtype="pdf", imgoutput="join",
      # Necessary for export
      exportpath=getwd())

p1 <- plotQ(alignList,imgoutput="join",returnplot=T,exportplot=F,basesize=11)
grid.arrange(p1$plot[[1]])

tr1 <- tabulateQ(qlist=fsList)
