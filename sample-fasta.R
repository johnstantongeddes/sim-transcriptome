## Script to select random fasta sequences from starting file

##-------------Load libraries--------------------------------
library(Biostrings)


##-------------Get command-line arguments and read files------------
args <- commandArgs(trailingOnly = TRUE)
filein <- args[1]
X <- args[2]

## Read known fasta file
file <- readDNAStringSet(filein)
length(file)

# Sample X sequences
samp.file <- sample(file, X)
length(samp.file)
head(samp.file)

# Write to file
writeXStringSet(samp.file, file="known.fasta")

