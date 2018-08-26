# elise_R_test2.R
# This script will compare before/after values for all plates to statistically identify difference

library(ggplot2)
library(gridExtra)
library(lme4)

# Read data
maindir = "E:/elise/R_scripts"
setwd(maindir)

# Specify output dir
outdir = paste(maindir, "/outputs", sep="")

# Perform statistics
## Test 5: growth before VERSUS growth after
### Find file
file = "quizcat_div_overall.csv"

### Perform test
quizData = read.csv(file)
stat = wilcox.test(quizData$div_before, quizData$div_after, paired=TRUE, alternative="greater")
#mdl = lm(count_category_before ~ count_category_after, data = quizData)

### Compute rough effect size
lessthan = 0
for (n in 1:nrow(quizData)){
  if (quizData$div_before[n] > quizData$div_after[n]){
    lessthan = lessthan + 1
  }
}
proportion = lessthan / nrow(quizData)
tmpTable = data.frame(stat = unlist(stat), stringsAsFactors = FALSE)
tmpTable[nrow(tmpTable) + 1,] = proportion
rownames(tmpTable)[7] = "proportion_less_than"

write.csv(tmpTable, file = paste(outdir,"/", "TEST5_div_beforeVSafter_wilcoxon.csv", sep=""))
          