# Script to calculate Transcripts per million (TPM) from gene expression count data
# From Wagner et al. 2012 "Measurement of mRNA abundance using RNA-seq data: RPKM measure
#  is inconsistent among samples" Theory Biosci. 131:281-285

library(plyr)
library(Biostrings)
library(stringr)

##-------------Get command-line arguments and read files------------
args <- commandArgs(trailingOnly = TRUE)
counts <- args[1]
rlsim <- args[2]

# read counts file
data <- read.table(counts)
# rename columns
colnames(data) <- c("locus","start","Length","count")
head(data)

cat("Done reading files", '\n')


##-------------Calculate TPM---------------------------------------

# TPM = (Rg * 10^6) / (Tn * Lg)
# where
# Rg: number of reads mapped to a particular transcript g = count
# Tn = sum of all length normalized transcript counts
# Lg: length of transcript g (kilobases)

# Calculate Tn
(Tn <- sum(ddply(data, .(locus), summarize, Tn = count/Length)[2]))

# Calculate TPM for each gene
tpm <- ddply(data, .(locus), summarize, tpm = (count*1e6)/(Tn*Length))
head(tpm)


##----------Compare to known expression levels---------------------

## Read in known expression levels
# read fasta file
known <- readDNAStringSet(rlsim)
# extract names
known.names <- names(known)
# strip expression level after $
kn.split <- str_split_fixed(known.names, '\\$', 2)
kn <- as.data.frame(kn.split)
colnames(kn) <- c("locus", "expression")
kn$expression <- as.numeric(as.character(kn$expression))
head(kn)
str(kn)

## Merge transcripts with both known and mapped expression
mdf <- merge(kn, tpm, by = "locus")
dim(mdf)
head(mdf)

## Correlation of known and assembled expression

cor.test(mdf$expression, mdf$tpm)

png("plot-known-vs-TPM-expression.png")
plot(mdf$expression, mdf$tpm, xlab="Known gene expression", ylab="Mapped gene expression")
dev.off()
