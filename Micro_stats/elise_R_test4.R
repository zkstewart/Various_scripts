# elise_R_stats.R
# This script will compute various statistical tests and produce figures for the microbiology study

library(ggplot2)
library(gridExtra)
library(FSA)

# Read data
maindir = "E:/elise/R_scripts"
setwd(maindir)

# Specify output dir
outdir = paste(maindir, "/outputs", sep="")

# Perform statistics
## Test 4: growth before / growth after for quiz responses
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
  else if (grepl('quizcat_div_before', val)){
    quizcat_files[[i]] = val
    i = i + 1
  }
}

### Loop through files and generate bar plots
quizcat_plots_test1 = list()
i = 1
for (file in quizcat_files){
  quizData = read.csv(file)
  plot = ggplot(quizData) + geom_bar(aes(x = factor(quiz_response), fill = factor(species_num)), position = "fill") + xlab("Quiz response number") + ylab("Proportion of answers") + ggtitle(paste("Quiz question #", toString(i))) + guides(fill=guide_legend(title="Species number"))
  quizcat_plots_test1[[i]] = plot
  i = i + 1
}
g=grid.arrange(grobs = quizcat_plots_test1[1:length(quizcat_files)])
ggsave(file=paste(outdir,"/", 'TEST4_quizcat_div_before.pdf', sep=""), g, width=25, height=10, units="in")

### Loop through files and perform statistics
quizcat_stats_test4 = list()
i = 1
for (file in quizcat_files){
  quizData = read.csv(file)
  fit = kruskal.test(species_num ~ quiz_response, data=quizData)
  fitaov = aov(species_num ~ quiz_response, data=quizData)
  quizcat_stats_test4[[i]] = fit
  i = i + 1
}
i = 1
for (stat in quizcat_stats_test4) {
  tmpTable = as.data.frame(stat[1:5])
  tmpTable$quiz_question = toString(i)
  if (i == 1) {
    mainDF = tmpTable
  }
  if (i != 1) {
    mainDF <- rbind(tmpTable[1, ], mainDF)
  }
  i = i + 1
}
### Not quite significance for quiz question 1! post hoc test time
quizData = read.csv(quizcat_files[[1]])
q1PH = dunnTest(species_num ~ quiz_response, data=quizData)
result = as.data.frame(q1PH$res)
result$a = 'NA'
result$b = 'NA'
tmpColNames = data.frame(as.list(colnames(result)))
### Do stupid stuff to make R not be terrible
mainDF <- data.frame(lapply(mainDF, as.character), stringsAsFactors=FALSE)
colnames(tmpColNames) = colnames(mainDF)
colnames(result) = colnames(mainDF)
mainDF <- rbind(colnames(mainDF), mainDF)
mainDF <- rbind(result, mainDF)
mainDF <- rbind(tmpColNames, mainDF)
### Write to file
write.table(mainDF, sep=",", file = paste(outdir,"/", "TEST4_quizcat_div_before_stats.csv", sep=""), row.names = FALSE, col.names = FALSE)

## Test 4: After
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
  else if (grepl('quizcat_div_after', val)){
    quizcat_files[[i]] = val
    i = i + 1
  }
}

### Loop through files and generate bar plots
quizcat_plots_test1 = list()
i = 1
for (file in quizcat_files){
  quizData = read.csv(file)
  plot = ggplot(quizData) + geom_bar(aes(x = factor(quiz_response), fill = factor(species_num)), position = "fill") + xlab("Quiz response number") + ylab("Proportion of answers") + ggtitle(paste("Quiz question #", toString(i))) + guides(fill=guide_legend(title="Species number"))
  quizcat_plots_test1[[i]] = plot
  i = i + 1
}
g=grid.arrange(grobs = quizcat_plots_test1[1:length(quizcat_files)])
ggsave(file=paste(outdir,"/", 'TEST4_quizcat_div_after.pdf',sep=""), g, width=25, height=10, units="in")

### Loop through files and perform statistics
quizcat_stats_test4 = list()
i = 1
for (file in quizcat_files){
  quizData = read.csv(file)
  fit = kruskal.test(species_num ~ quiz_response, data=quizData)
  quizcat_stats_test4[[i]] = fit
  i = i + 1
}
i = 1
for (stat in quizcat_stats_test4) {
  tmpTable = as.data.frame(stat[1:5])
  tmpTable$quiz_question = toString(i)
  if (i == 1) {
    mainDF = tmpTable
  }
  if (i != 1) {
    mainDF <- rbind(tmpTable[1, ], mainDF)
  }
  i = i + 1
}
write.csv(mainDF, file = paste(outdir,"/", "TEST4_quizcat_div_after_stats.csv", sep=""))
