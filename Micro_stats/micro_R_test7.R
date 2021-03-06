# elise_R_stats.R
# This script will compute various statistical tests and produce figures for the microbiology study

library(ggplot2)
library(gridExtra)
library(FSA)
library(MASS)
library(rcompanion)

# Read data
maindir = "E:/elise/R_scripts"
setwd(maindir)

# Specify output dir
outdir = paste(maindir, "/outputs", sep="")

# Perform statistics
## Test 7: growth before / growth after for quiz responses
### Find files
currdir = list.files()
quizcat_files_NB = list()
quizcat_files_HBA = list()
quizcat_files_MSA = list()
quizcat_files_MCA = list()
#### Loop NB
i = 1
for (val in currdir){
  if(grepl('.pdf', val)){
    donothing=1
  } else if (grepl('_stats', val)){
    donothing=1
  }
  else if (grepl('quizcat_media_before_NB', val)){
    quizcat_files_NB[[i]] = val
    i = i + 1
  }
}
#### Loop HBA
i = 1
for (val in currdir){
  if(grepl('.pdf', val)){
    donothing=1
  } else if (grepl('_stats', val)){
    donothing=1
  }
  else if (grepl('quizcat_media_before_HBA', val)){
    quizcat_files_HBA[[i]] = val
    i = i + 1
  }
}
#### Loop MSA
i = 1
for (val in currdir){
  if(grepl('.pdf', val)){
    donothing=1
  } else if (grepl('_stats', val)){
    donothing=1
  }
  else if (grepl('quizcat_media_before_MSA', val)){
    quizcat_files_MSA[[i]] = val
    i = i + 1
  }
}
#### Loop MCA
i = 1
for (val in currdir){
  if(grepl('.pdf', val)){
    donothing=1
  } else if (grepl('_stats', val)){
    donothing=1
  }
  else if (grepl('quizcat_media_before_MCA', val)){
    quizcat_files_MCA[[i]] = val
    i = i + 1
  }
}


### Loop through files and generate bar plots
quizcat_plots_test7NB = list()
quizcat_plots_test7HBA = list()
quizcat_plots_test7MSA = list()
quizcat_plots_test7MCA = list()
#### NB LOOP
i = 1
for (file in quizcat_files_NB){
  quizData = read.csv(file)
  plot = ggplot(quizData) + geom_bar(aes(x = factor(quiz_response), fill = factor(reaction_cat)), position = "fill") + xlab("Quiz response number") + ylab("Proportion of answers") + ggtitle(paste("Quiz question #", toString(i))) + guides(fill=guide_legend(title="Number of reaction type"))
  plot = plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  quizcat_plots_test7NB[[i]] = plot
  i = i + 1
}
g=grid.arrange(grobs = quizcat_plots_test7NB[1:length(quizcat_files_NB)], ncol = 2)
ggsave(file=paste(outdir,"/", 'TEST7_quizcat_NB_before.pdf', sep=""), g, width=15, height=10, units="in")

#### HBA LOOP
i = 1
for (file in quizcat_files_HBA){
  quizData = read.csv(file)
  plot = ggplot(quizData) + geom_bar(aes(x = factor(quiz_response), fill = factor(reaction_cat)), position = "fill") + xlab("Quiz response number") + ylab("Proportion of answers") + ggtitle(paste("Quiz question #", toString(i))) + guides(fill=guide_legend(title="Number of reaction type"))
  plot = plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  quizcat_plots_test7HBA[[i]] = plot
  i = i + 1
}
g=grid.arrange(grobs = quizcat_plots_test7HBA[1:length(quizcat_files_HBA)], ncol = 2)
ggsave(file=paste(outdir,"/", 'TEST7_quizcat_HBA_before.pdf', sep=""), g, width=15, height=10, units="in")

#### MSA LOOP
i = 1
for (file in quizcat_files_MSA){
  quizData = read.csv(file)
  plot = ggplot(quizData) + geom_bar(aes(x = factor(quiz_response), fill = factor(reaction_cat)), position = "fill") + xlab("Quiz response number") + ylab("Proportion of answers") + ggtitle(paste("Quiz question #", toString(i))) + guides(fill=guide_legend(title="Number of reaction type"))
  plot = plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  quizcat_plots_test7MSA[[i]] = plot
  i = i + 1
}
g=grid.arrange(grobs = quizcat_plots_test7MSA[1:length(quizcat_files_MSA)], ncol = 2)
ggsave(file=paste(outdir,"/", 'TEST7_quizcat_MSA_before.pdf', sep=""), g, width=15, height=10, units="in")

#### MCA LOOP
i = 1
for (file in quizcat_files_MCA){
  quizData = read.csv(file)
  plot = ggplot(quizData) + geom_bar(aes(x = factor(quiz_response), fill = factor(reaction_cat)), position = "fill") + xlab("Quiz response number") + ylab("Proportion of answers") + ggtitle(paste("Quiz question #", toString(i))) + guides(fill=guide_legend(title="Number of reaction type"))
  plot = plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  quizcat_plots_test7MCA[[i]] = plot
  i = i + 1
}
g=grid.arrange(grobs = quizcat_plots_test7MCA[1:length(quizcat_files_MCA)], ncol = 2)
ggsave(file=paste(outdir,"/", 'TEST7_quizcat_MCA_before.pdf', sep=""), g, width=15, height=10, units="in")


### Loop through files and perform statistics
quizcat_stats_test7NB = list()
quizcat_stats_test7HBA = list()
quizcat_stats_test7MSA = list()
quizcat_stats_test7MCA = list()
#### NB LOOP
i = 1
for (file in quizcat_files_NB){
  quizData = read.csv(file)
  tbl = table(quizData$quiz_response, quizData$reaction_cat)
  #fit = fisher.test(tbl,workspace=2e+07)
  fit = chisq.test(tbl)
  #fit = kruskal.test(reaction_cat ~ quiz_response, data=quizData)
  quizcat_stats_test7NB[[i]] = fit
  i = i + 1
}
i = 1
for (stat in quizcat_stats_test7NB) {
  tmpTable = as.data.frame(stat[1:4])
  tmpTable$quiz_question = toString(i)
  if (i == 1) {
    mainDF = tmpTable
  }
  if (i != 1) {
    mainDF <- rbind(tmpTable[1, ], mainDF)
  }
  i = i + 1
} #8
          # Significance for quiz question 8! post hoc test time
quizData = read.csv(quizcat_files_NB[[8]])
tbl = table(quizData$quiz_response, quizData$reaction_cat)
q1PH = pairwiseNominalIndependence(tbl,fisher = TRUE,gtest = FALSE,chisq = FALSE,digits = 3)
result = as.data.frame(q1PH)
result$quiz_q_8 = 'NA'
result$b = 'NA'
result = rbind(colnames(result), result)
colnames(result) = c("p.value","alternative","method","data.name","quiz_question")
tmpColNames = data.frame(as.list(colnames(result)))
### Do stupid stuff to make R not be terrible
mainDF <- data.frame(lapply(mainDF, as.character), stringsAsFactors=FALSE)
colnames(tmpColNames) = colnames(mainDF)
colnames(result) = colnames(mainDF)
mainDF <- rbind(colnames(mainDF), mainDF)
mainDF <- rbind(result, mainDF)

write.table(mainDF, sep=",", file = paste(outdir,"/", "TEST7_quizcat_NB_before_stats.csv", sep=""), row.names = FALSE, col.names = FALSE)


#### HBA LOOP
i = 1
for (file in quizcat_files_HBA){
  quizData = read.csv(file)
  tbl = table(quizData$quiz_response, quizData$reaction_cat)
  #fit = fisher.test(tbl,workspace=2e+07,hybrid = TRUE)
  fit = chisq.test(tbl)
  #fit = kruskal.test(reaction_cat ~ quiz_response, data=quizData)
  quizcat_stats_test7HBA[[i]] = fit
  i = i + 1
}
i = 1
for (stat in quizcat_stats_test7HBA) {
  tmpTable = as.data.frame(stat[1:4])
  tmpTable$quiz_question = toString(i)
  if (i == 1) {
    mainDF = tmpTable
  }
  if (i != 1) {
    mainDF <- rbind(tmpTable[1, ], mainDF)
  }
  i = i + 1
} #1,3,5,6,8
# Significance for quiz question 1! post hoc test time
quizData = read.csv(quizcat_files_HBA[[1]])
tbl = table(quizData$quiz_response, quizData$reaction_cat)
q1PH = pairwiseNominalIndependence(tbl,fisher = TRUE,gtest = FALSE,chisq = FALSE,digits = 3)
result = as.data.frame(q1PH)
result$quiz_q_1 = 'NA'
result$b = 'NA'
result = rbind(colnames(result), result)
colnames(result) = c("p.value","alternative","method","data.name","quiz_question")
tmpColNames = data.frame(as.list(colnames(result)))
### Do stupid stuff to make R not be terrible
mainDF <- data.frame(lapply(mainDF, as.character), stringsAsFactors=FALSE)
colnames(tmpColNames) = colnames(mainDF)
colnames(result) = colnames(mainDF)
mainDF <- rbind(colnames(mainDF), mainDF)
mainDF <- rbind(result, mainDF)
# Significance for quiz question 3! post hoc test time
quizData = read.csv(quizcat_files_HBA[[3]])
tbl = table(quizData$quiz_response, quizData$reaction_cat)
q1PH = pairwiseNominalIndependence(tbl,fisher = TRUE,gtest = FALSE,chisq = FALSE,digits = 3)
result = as.data.frame(q1PH)
result$quiz_q_3 = 'NA'
result$b = 'NA'
result = rbind(colnames(result), result)
colnames(result) = c("p.value","alternative","method","data.name","quiz_question")
tmpColNames = data.frame(as.list(colnames(result)))
### Do stupid stuff to make R not be terrible
mainDF <- data.frame(lapply(mainDF, as.character), stringsAsFactors=FALSE)
colnames(tmpColNames) = colnames(mainDF)
colnames(result) = colnames(mainDF)
#mainDF <- rbind(colnames(mainDF), mainDF)
mainDF <- rbind(result, mainDF)
            # Significance for quiz question 5! post hoc test time
quizData = read.csv(quizcat_files_HBA[[5]])
tbl = table(quizData$quiz_response, quizData$reaction_cat)
q1PH = pairwiseNominalIndependence(tbl,fisher = TRUE,gtest = FALSE,chisq = FALSE,digits = 3)
result = as.data.frame(q1PH)
result$quiz_q_5 = 'NA'
result$b = 'NA'
result = rbind(colnames(result), result)
colnames(result) = c("p.value","alternative","method","data.name","quiz_question")
tmpColNames = data.frame(as.list(colnames(result)))
### Do stupid stuff to make R not be terrible
mainDF <- data.frame(lapply(mainDF, as.character), stringsAsFactors=FALSE)
colnames(tmpColNames) = colnames(mainDF)
colnames(result) = colnames(mainDF)
#mainDF <- rbind(colnames(mainDF), mainDF)
mainDF <- rbind(result, mainDF)
          # Significance for quiz question 6! post hoc test time
quizData = read.csv(quizcat_files_HBA[[6]])
tbl = table(quizData$quiz_response, quizData$reaction_cat)
q1PH = pairwiseNominalIndependence(tbl,fisher = TRUE,gtest = FALSE,chisq = FALSE,digits = 3)
result = as.data.frame(q1PH)
result$quiz_q_6 = 'NA'
result$b = 'NA'
result = rbind(colnames(result), result)
colnames(result) = c("p.value","alternative","method","data.name","quiz_question")
tmpColNames = data.frame(as.list(colnames(result)))
### Do stupid stuff to make R not be terrible
mainDF <- data.frame(lapply(mainDF, as.character), stringsAsFactors=FALSE)
colnames(tmpColNames) = colnames(mainDF)
colnames(result) = colnames(mainDF)
#mainDF <- rbind(colnames(mainDF), mainDF)
mainDF <- rbind(result, mainDF)
        # Significance for quiz question 8! post hoc test time
quizData = read.csv(quizcat_files_HBA[[8]])
tbl = table(quizData$quiz_response, quizData$reaction_cat)
q1PH = pairwiseNominalIndependence(tbl,fisher = TRUE,gtest = FALSE,chisq = FALSE,digits = 3)
result = as.data.frame(q1PH)
result$quiz_q_8 = 'NA'
result$b = 'NA'
result = rbind(colnames(result), result)
colnames(result) = c("p.value","alternative","method","data.name","quiz_question")
tmpColNames = data.frame(as.list(colnames(result)))
### Do stupid stuff to make R not be terrible
mainDF <- data.frame(lapply(mainDF, as.character), stringsAsFactors=FALSE)
colnames(tmpColNames) = colnames(mainDF)
colnames(result) = colnames(mainDF)
#mainDF <- rbind(colnames(mainDF), mainDF)
mainDF <- rbind(result, mainDF)
### Write to file
write.table(mainDF, sep=",", file = paste(outdir,"/", "TEST7_quizcat_HBA_before_stats.csv", sep=""), row.names = FALSE, col.names = FALSE)

#### MSA LOOP
i = 1
for (file in quizcat_files_MSA){
  quizData = read.csv(file)
  tbl = table(quizData$quiz_response, quizData$reaction_cat)
  #fit = fisher.test(tbl,workspace=2e+07)
  fit = chisq.test(tbl)
  #fit = kruskal.test(reaction_cat ~ quiz_response, data=quizData)
  quizcat_stats_test7MSA[[i]] = fit
  i = i + 1
}
i = 1
for (stat in quizcat_stats_test7MSA) {
  tmpTable = as.data.frame(stat[1:4])
  tmpTable$quiz_question = toString(i)
  if (i == 1) {
    mainDF = tmpTable
  }
  if (i != 1) {
    mainDF <- rbind(tmpTable[1, ], mainDF)
  }
  i = i + 1
} #1,5
# Significance for quiz question 1! post hoc test time
quizData = read.csv(quizcat_files_MSA[[1]])
tbl = table(quizData$quiz_response, quizData$reaction_cat)
q1PH = pairwiseNominalIndependence(tbl,fisher = TRUE,gtest = FALSE,chisq = FALSE,digits = 3)
result = as.data.frame(q1PH)
result$quiz_q_1 = 'NA'
result$b = 'NA'
result = rbind(colnames(result), result)
colnames(result) = c("p.value","alternative","method","data.name","quiz_question")
#tmpColNames = data.frame(as.list(colnames(result)))
### Do stupid stuff to make R not be terrible
mainDF <- data.frame(lapply(mainDF, as.character), stringsAsFactors=FALSE)
#colnames(tmpColNames) = colnames(mainDF)
colnames(result) = colnames(mainDF)
mainDF <- rbind(colnames(mainDF), mainDF)
mainDF <- rbind(result, mainDF)
# Significance for quiz question 5! post hoc test time
quizData = read.csv(quizcat_files_MSA[[5]])
tbl = table(quizData$quiz_response, quizData$reaction_cat)
q1PH = pairwiseNominalIndependence(tbl,fisher = TRUE,gtest = FALSE,chisq = FALSE,digits = 3)
result = as.data.frame(q1PH)
result$quiz_q_5 = 'NA'
result$b = 'NA'
result = rbind(colnames(result), result)
colnames(result) = c("p.value","alternative","method","data.name","quiz_question")
tmpColNames = data.frame(as.list(colnames(result)))
### Do stupid stuff to make R not be terrible
mainDF <- data.frame(lapply(mainDF, as.character), stringsAsFactors=FALSE)
colnames(tmpColNames) = colnames(mainDF)
colnames(result) = colnames(mainDF)
#mainDF <- rbind(colnames(mainDF), mainDF)
mainDF <- rbind(result, mainDF)
### Write to file
write.table(mainDF, sep=",", file = paste(outdir,"/", "TEST7_quizcat_MSA_before_stats.csv", sep=""), row.names = FALSE, col.names = FALSE)


#### MCA LOOP
i = 1
for (file in quizcat_files_MCA){
  quizData = read.csv(file)
  tbl = table(quizData$quiz_response, quizData$reaction_cat)
  #fit = fisher.test(tbl,workspace=2e+07)
  fit = chisq.test(tbl)
  #fit = kruskal.test(reaction_cat ~ quiz_response, data=quizData)
  quizcat_stats_test7MCA[[i]] = fit
  i = i + 1
}
i = 1
for (stat in quizcat_stats_test7MCA) {
  tmpTable = as.data.frame(stat[1:4])
  tmpTable$quiz_question = toString(i)
  if (i == 1) {
    mainDF = tmpTable
  }
  if (i != 1) {
    mainDF <- rbind(tmpTable[1, ], mainDF)
  }
  i = i + 1
}
write.table(mainDF, sep=",", file = paste(outdir,"/", "TEST7_quizcat_MCA_before_stats.csv", sep=""), row.names = FALSE, col.names = FALSE)

## Test 7: growth after for quiz responses
### Find files
currdir = list.files()
quizcat_files_NB = list()
quizcat_files_HBA = list()
quizcat_files_MSA = list()
quizcat_files_MCA = list()
#### Loop NB
i = 1
for (val in currdir){
  if(grepl('.pdf', val)){
    donothing=1
  } else if (grepl('_stats', val)){
    donothing=1
  }
  else if (grepl('quizcat_media_after_NB', val)){
    quizcat_files_NB[[i]] = val
    i = i + 1
  }
}
#### Loop HBA
i = 1
for (val in currdir){
  if(grepl('.pdf', val)){
    donothing=1
  } else if (grepl('_stats', val)){
    donothing=1
  }
  else if (grepl('quizcat_media_after_HBA', val)){
    quizcat_files_HBA[[i]] = val
    i = i + 1
  }
}
#### Loop MSA
i = 1
for (val in currdir){
  if(grepl('.pdf', val)){
    donothing=1
  } else if (grepl('_stats', val)){
    donothing=1
  }
  else if (grepl('quizcat_media_after_MSA', val)){
    quizcat_files_MSA[[i]] = val
    i = i + 1
  }
}
#### Loop MCA
i = 1
for (val in currdir){
  if(grepl('.pdf', val)){
    donothing=1
  } else if (grepl('_stats', val)){
    donothing=1
  }
  else if (grepl('quizcat_media_after_MCA', val)){
    quizcat_files_MCA[[i]] = val
    i = i + 1
  }
}

### Loop through files and generate bar plots
quizcat_plots_test7NB = list()
quizcat_plots_test7HBA = list()
quizcat_plots_test7MSA = list()
quizcat_plots_test7MCA = list()
#### NB LOOP
i = 1
for (file in quizcat_files_NB){
  quizData = read.csv(file)
  plot = ggplot(quizData) + geom_bar(aes(x = factor(quiz_response), fill = factor(reaction_cat)), position = "fill") + xlab("Quiz response number") + ylab("Proportion of answers") + ggtitle(paste("Quiz question #", toString(i))) + guides(fill=guide_legend(title="Number of reaction type"))
  plot = plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  quizcat_plots_test7NB[[i]] = plot
  i = i + 1
}
g=grid.arrange(grobs = quizcat_plots_test7NB[1:length(quizcat_files_NB)], ncol = 2)
ggsave(file=paste(outdir,"/", 'TEST7_quizcat_NB_after.pdf', sep=""), g, width=15, height=10, units="in")

#### HBA LOOP
i = 1
for (file in quizcat_files_HBA){
  quizData = read.csv(file)
  plot = ggplot(quizData) + geom_bar(aes(x = factor(quiz_response), fill = factor(reaction_cat)), position = "fill") + xlab("Quiz response number") + ylab("Proportion of answers") + ggtitle(paste("Quiz question #", toString(i))) + guides(fill=guide_legend(title="Number of reaction type"))
  plot = plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  quizcat_plots_test7HBA[[i]] = plot
  i = i + 1
}
g=grid.arrange(grobs = quizcat_plots_test7HBA[1:length(quizcat_files_HBA)], ncol = 2)
ggsave(file=paste(outdir,"/", 'TEST7_quizcat_HBA_after.pdf', sep=""), g, width=15, height=10, units="in")

#### MSA LOOP
i = 1
for (file in quizcat_files_MSA){
  quizData = read.csv(file)
  plot = ggplot(quizData) + geom_bar(aes(x = factor(quiz_response), fill = factor(reaction_cat)), position = "fill") + xlab("Quiz response number") + ylab("Proportion of answers") + ggtitle(paste("Quiz question #", toString(i))) + guides(fill=guide_legend(title="Number of reaction type"))
  plot = plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  quizcat_plots_test7MSA[[i]] = plot
  i = i + 1
}
g=grid.arrange(grobs = quizcat_plots_test7MSA[1:length(quizcat_files_MSA)], ncol = 2)
ggsave(file=paste(outdir,"/", 'TEST7_quizcat_MSA_after.pdf', sep=""), g, width=15, height=10, units="in")

#### MCA LOOP
i = 1
for (file in quizcat_files_MCA){
  quizData = read.csv(file)
  plot = ggplot(quizData) + geom_bar(aes(x = factor(quiz_response), fill = factor(reaction_cat)), position = "fill") + xlab("Quiz response number") + ylab("Proportion of answers") + ggtitle(paste("Quiz question #", toString(i))) + guides(fill=guide_legend(title="Number of reaction type"))
  plot = plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  quizcat_plots_test7MCA[[i]] = plot
  i = i + 1
}
g=grid.arrange(grobs = quizcat_plots_test7MCA[1:length(quizcat_files_MCA)], ncol = 2)
ggsave(file=paste(outdir,"/", 'TEST7_quizcat_MCA_after.pdf', sep=""), g, width=15, height=10, units="in")


### Loop through files and perform statistics
quizcat_stats_test7NB = list()
quizcat_stats_test7HBA = list()
quizcat_stats_test7MSA = list()
quizcat_stats_test7MCA = list()
#### NB LOOP
i = 1
for (file in quizcat_files_NB){
  quizData = read.csv(file)
  tbl = table(quizData$quiz_response, quizData$reaction_cat)
  #fit = fisher.test(tbl,workspace=2e+07)
  fit = chisq.test(tbl)
  #fit = kruskal.test(reaction_cat ~ quiz_response, data=quizData)
  quizcat_stats_test7NB[[i]] = fit
  i = i + 1
}
i = 1
for (stat in quizcat_stats_test7NB) {
  tmpTable = as.data.frame(stat[1:4])
  tmpTable$quiz_question = toString(i)
  if (i == 1) {
    mainDF = tmpTable
  }
  if (i != 1) {
    mainDF <- rbind(tmpTable[1, ], mainDF)
  }
  i = i + 1
} #4
# Significance for quiz question 4! post hoc test time
quizData = read.csv(quizcat_files_NB[[4]])
tbl = table(quizData$quiz_response, quizData$reaction_cat)
q1PH = pairwiseNominalIndependence(tbl,fisher = TRUE,gtest = FALSE,chisq = FALSE,digits = 3)
result = as.data.frame(q1PH)
result$quiz_q_4 = 'NA'
result$b = 'NA'
result = rbind(colnames(result), result)
colnames(result) = c("p.value","alternative","method","data.name","quiz_question")
tmpColNames = data.frame(as.list(colnames(result)))
### Do stupid stuff to make R not be terrible
mainDF <- data.frame(lapply(mainDF, as.character), stringsAsFactors=FALSE)
colnames(tmpColNames) = colnames(mainDF)
colnames(result) = colnames(mainDF)
mainDF <- rbind(colnames(mainDF), mainDF)
mainDF <- rbind(result, mainDF)
### Write to file
write.table(mainDF, sep=",", file = paste(outdir,"/", "TEST7_quizcat_NB_after_stats.csv", sep=""), row.names = FALSE, col.names = FALSE)

#### HBA LOOP
i = 1
for (file in quizcat_files_HBA){
  quizData = read.csv(file)
  tbl = table(quizData$quiz_response, quizData$reaction_cat)
  #fit = fisher.test(tbl,workspace=2e+07)
  fit = chisq.test(tbl)
  #fit = kruskal.test(reaction_cat ~ quiz_response, data=quizData)
  quizcat_stats_test7HBA[[i]] = fit
  i = i + 1
}
i = 1
for (stat in quizcat_stats_test7HBA) {
  tmpTable = as.data.frame(stat[1:4])
  tmpTable$quiz_question = toString(i)
  if (i == 1) {
    mainDF = tmpTable
  }
  if (i != 1) {
    mainDF <- rbind(tmpTable[1, ], mainDF)
  }
  i = i + 1
} #6
# Significance for quiz question 6! post hoc test time
quizData = read.csv(quizcat_files_HBA[[6]])
tbl = table(quizData$quiz_response, quizData$reaction_cat)
q1PH = pairwiseNominalIndependence(tbl,fisher = TRUE,gtest = FALSE,chisq = FALSE,digits = 3)
result = as.data.frame(q1PH)
result$quiz_q_6 = 'NA'
result$b = 'NA'
result = rbind(colnames(result), result)
colnames(result) = c("p.value","alternative","method","data.name","quiz_question")
tmpColNames = data.frame(as.list(colnames(result)))
### Do stupid stuff to make R not be terrible
mainDF <- data.frame(lapply(mainDF, as.character), stringsAsFactors=FALSE)
colnames(tmpColNames) = colnames(mainDF)
colnames(result) = colnames(mainDF)
mainDF <- rbind(colnames(mainDF), mainDF)
mainDF <- rbind(result, mainDF)
### Write to file
write.table(mainDF, sep=",", file = paste(outdir,"/", "TEST7_quizcat_HBA_after_stats.csv", sep=""), row.names = FALSE, col.names = FALSE)

#### MSA LOOP
i = 1
for (file in quizcat_files_MSA){
  quizData = read.csv(file)
  tbl = table(quizData$quiz_response, quizData$reaction_cat)
  #fit = fisher.test(tbl,workspace=2e+07)
  fit = chisq.test(tbl)
  #fit = kruskal.test(reaction_cat ~ quiz_response, data=quizData)
  quizcat_stats_test7MSA[[i]] = fit
  i = i + 1
}
i = 1
for (stat in quizcat_stats_test7MSA) {
  tmpTable = as.data.frame(stat[1:4])
  tmpTable$quiz_question = toString(i)
  if (i == 1) {
    mainDF = tmpTable
  }
  if (i != 1) {
    mainDF <- rbind(tmpTable[1, ], mainDF)
  }
  i = i + 1
}
### Write to file
write.table(mainDF, sep=",", file = paste(outdir,"/", "TEST7_quizcat_MSA_after_stats.csv", sep=""), row.names = FALSE, col.names = FALSE)

#### MCA LOOP
i = 1
for (file in quizcat_files_MCA){
  quizData = read.csv(file)
  tbl = table(quizData$quiz_response, quizData$reaction_cat)
  #fit = fisher.test(tbl,workspace=2e+07)
  fit = chisq.test(tbl)
  #fit = kruskal.test(reaction_cat ~ quiz_response, data=quizData)
  quizcat_stats_test7MCA[[i]] = fit
  i = i + 1
}
i = 1
for (stat in quizcat_stats_test7MCA) {
  tmpTable = as.data.frame(stat[1:4])
  tmpTable$quiz_question = toString(i)
  if (i == 1) {
    mainDF = tmpTable
  }
  if (i != 1) {
    mainDF <- rbind(tmpTable[1, ], mainDF)
  }
  i = i + 1
} #4
# Significance for quiz question 4! post hoc test time
quizData = read.csv(quizcat_files_MCA[[4]])
tbl = table(quizData$quiz_response, quizData$reaction_cat)
q1PH = pairwiseNominalIndependence(tbl,fisher = TRUE,gtest = FALSE,chisq = FALSE,digits = 3)
result = as.data.frame(q1PH)
result$quiz_q_4 = 'NA'
result$b = 'NA'
result = rbind(colnames(result), result)
colnames(result) = c("p.value","alternative","method","data.name","quiz_question")
tmpColNames = data.frame(as.list(colnames(result)))
### Do stupid stuff to make R not be terrible
mainDF <- data.frame(lapply(mainDF, as.character), stringsAsFactors=FALSE)
colnames(tmpColNames) = colnames(mainDF)
colnames(result) = colnames(mainDF)
mainDF <- rbind(colnames(mainDF), mainDF)
mainDF <- rbind(result, mainDF)
### Write to file
write.table(mainDF, sep=",", file = paste(outdir,"/", "TEST7_quizcat_MCA_after_stats.csv", sep=""), row.names = FALSE, col.names = FALSE)
