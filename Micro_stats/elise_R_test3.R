# elise_R_stats.R
# This script will compute various statistical tests and produce figures for the microbiology study

library(ggplot2)
library(gridExtra)
library(lme4)
#library(devtools)
#install_github("dgrtwo/broom")
library(broom)

# Read data
maindir = "E:/elise/R_scripts"
setwd(maindir)

# Specify output dir
outdir = paste(maindir, "/outputs", sep="")

# Perform statistics
## Test 3: growth before VERSUS growth after for quiz responses
### Find files
currdir = list.files()
i = 1
quizcat_files = list()
for (val in currdir){
  if(grepl('.pdf', val)){
    donothing=1
  } else if (grepl('_stats', val)){
    donothing=1
  }
  else if (grepl('quizcat_totals_compare', val)){
    quizcat_files[[i]] = val
    i = i + 1
  }
}

### Loop through files and perform statistics
quizcat_wcx_test2 = list()
quizcat_lm_test2 = list()

i = 1
for (file in quizcat_files){
  quizData = read.csv(file)
  #fit = wilcox.test(quizData$count_category_before, quizData$count_category_after, paired=TRUE, alternative="greater")
  mdl = lm(count_category ~ treatment*quiz_response, data = quizData)
  #quizcat_wcx_test2[[i]] = fit
  quizcat_lm_test2[[i]] = mdl
  i = i + 1
}

## Linear model output - should be done with combined responses!
i = 1
for (stat in quizcat_lm_test2) {
  #tmpTable = as.data.frame(stat[1:7])
  tmpTable = tidy(stat)
  tmpTable[nrow(tmpTable) + 1,] = quizcat_files[[i]]
  tmpTable[nrow(tmpTable) + 1,] = ''
  if (i == 1) {
    mainDF = tmpTable
  }
  if (i != 1) {
    mainDF <- rbind(tmpTable, mainDF)
  }
  i = i + 1
}
write.csv(mainDF, file = paste(outdir,"/", "TEST3_quizcat_totals_compare_linear_stats.csv", sep=""), row.names = FALSE)

## Wilcoxon output - should be done with response separation!
#i = 1
#for (stat in quizcat_wcx_test2) {
#  #tmpTable = as.data.frame(stat[1:7])
#  tmpTable = data.frame(stat = unlist(stat), stringsAsFactors = FALSE)
#  tmpTable[nrow(tmpTable) + 1,] = quizcat_files[[i]]
#  rownames(tmpTable)[6] = "quiz_question"
#  if (i == 1) {
#    mainDF = tmpTable
#  }
#  if (i != 1) {
#    mainDF <- cbind(tmpTable, mainDF)
#  }
#  i = i + 1
#}
#outDF = t(mainDF)
#write.csv(outDF, file = paste(outdir,"/", "TEST3_quizcat_totals_compare_wilcoxon_stats.csv",sep=""))