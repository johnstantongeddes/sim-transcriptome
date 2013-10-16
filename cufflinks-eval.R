#################################################################
## Script to evaluate gene counts from cufflinks on known vs.
## assembled transcripts
##
## John Stanton-Geddes
## 2013-10-16
#################################################################


##-------------------Options-------------------------------------
library(stringr)
options(stringsAsFactors = FALSE)

##-------------Get command-line arguments and read files------------
args <- commandArgs(trailingOnly = TRUE)
data <- args[1]

# Read file
data <- read.table(data, header=TRUE, sep="\t")
head(data)


##---------Compare gene expression counts mapped to known transcripts with 'truth'---
split <- str_split_fixed(data$locus, ":", 2)
head(split)
split2 <- str_split_fixed(split[,1], '\\$', 2)
head(split2)

starting.expression <- as.numeric(split2[,2])
starting.expression

cor.test(starting.expression, data$FPKM)

pdf("corr-known-vs-assembled-expression-cufflinks.pdf")
plot(starting.expression, data$FPKM, xlab="Known expression", ylab="Assembled expression")
dev.off()




