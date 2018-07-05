# elise_R_stats.R
# This script will compute various statistical tests and produce figures for the microbiology study

library(ggplot2)
library(gridExtra)

# Read data
setwd("E:/elise/R_scripts")

# Perform statistics
## Test 1: growth before / growth after for quiz responses
### Find files
currdir = list.files()
i = 1
quizcat_files = list()
for (val in currdir){
  if(grepl('quizcat_totals_before', val)){
    quizcat_files[[i]] = val
    i = i + 1
  }
}

### Loop through files and generate bar plots
quizcat_plots_test1 = list()
i = 1
for (file in quizcat_files){
  quizData = read.csv(file)
  quizData$count_category <- factor(quizData$count_category, levels = c("0", "1-10", "11-100", "100-340", "TNTC"))
  plot = ggplot(quizData) + geom_bar(aes(x = factor(quiz_response), fill = factor(count_category)), position = "fill") + xlab("Quiz response number") + ylab("Proportion of answers") + ggtitle(paste("Quiz question #", toString(i+1))) + guides(fill=guide_legend(title="Size category"))
  quizcat_plots_test1[[i]] = plot
  i = i + 1
}
g=grid.arrange(grobs = quizcat_plots_test1[1:length(quizcat_files)])
ggsave(file='quizcat_totals_before.pdf', g, width=15, height=10, units="in")

### Loop through files and perform statistics
quizcat_stats_test1 = list()
i = 1
for (file in quizcat_files){
  quizData = read.csv(file)
  fit = kruskal.test(count_category ~ quiz_response, data=quizData)
  quizcat_stats_test1[[i]] = fit
  i = i + 1
}
i = 1
for (stat in quizcat_stats_test1) {
  tmpTable = as.data.frame(stat[1:5])
  tmpTable$quiz_question = toString(i+1)
  if (i == 1) {
    mainDF = tmpTable
  }
  if (i != 1) {
    mainDF <- rbind(tmpTable[1, ], mainDF)
  }
  i = i + 1
}
write.csv(mainDF, file = "quizcat_totals_before_stats.csv")

### Find files
currdir = list.files()
i = 1
quizcat_files = list()
for (val in currdir){
  if(grepl('quizcat_totals_after', val)){
    quizcat_files[[i]] = val
    i = i + 1
  }
}

### Loop through files and generate bar plots
quizcat_plots_test1 = list()
i = 1
for (file in quizcat_files){
  quizData = read.csv(file)
  quizData$count_category <- factor(quizData$count_category, levels = c("0", "1-10", "11-100", "100-340", "TNTC"))
  plot = ggplot(quizData) + geom_bar(aes(x = factor(quiz_response), fill = factor(count_category)), position = "fill") + xlab("Quiz response number") + ylab("Proportion of answers") + ggtitle(paste("Quiz question #", toString(i+1))) + guides(fill=guide_legend(title="Size category"))
  quizcat_plots_test1[[i]] = plot
  i = i + 1
}
g=grid.arrange(grobs = quizcat_plots_test1[1:length(quizcat_files)])
ggsave(file='quizcat_totals_after.pdf', g, width=15, height=10, units="in")

### Loop through files and perform statistics
quizcat_stats_test1 = list()
i = 1
for (file in quizcat_files){
  quizData = read.csv(file)
  fit = kruskal.test(count_category ~ quiz_response, data=quizData)
  quizcat_stats_test1[[i]] = fit
  i = i + 1
}
i = 1
for (stat in quizcat_stats_test1) {
  tmpTable = as.data.frame(stat[1:5])
  tmpTable$quiz_question = toString(i+1)
  if (i == 1) {
    mainDF = tmpTable
  }
  if (i != 1) {
    mainDF <- rbind(tmpTable[1, ], mainDF)
  }
  i = i + 1
}
write.csv(mainDF, file = "quizcat_totals_after_stats.csv")
