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
mean(width(file))

# Sample X sequences
samp.file <- sample(file, X)
length(samp.file)
head(samp.file)

# Write subset of fasta sequences to file
writeXStringSet(samp.file, file="known.fasta")

# Write mean length to file
cat(mean(width(samp.file)), '\n', file="mean-length-fasta.txt")
