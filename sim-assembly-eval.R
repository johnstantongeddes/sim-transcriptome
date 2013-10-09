#####################################################
## Script to evaluate simulated transcriptome assembly
##
## Takes in list of known genes and results of BLAST
## of known genes against transcriptome assembly
##
## Reports the starting number of genes, number of genes
## captured in assembly, number of genes missing in assembly
## and the average length of each of genes in these groups.
## Also produces histogram showing distribution of captured
## versus missing genes.
##
## Script expects two arguments:
##
##    Rscript sim-assembly-eval.R fasta blastout
##
## where fasta is the starting fasta file and
## blastout is the result of running blastn of the
## original fasta file on the assembled fasta specifying
## `-outfmt 6` to allow reading with read.blast() 
## 
## Author: John Stanton-Geddes
## 2013-10-09
#####################################################


##-------------Load libraries--------------------------------
#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
library(Biostrings)
library(stringr)
library(RFLPtools)
library(ggplot2)

##-------------Get command-line arguments and read files------------
args <- commandArgs(trailingOnly = TRUE)
fasta <- args[1]
blastout <- args[2]

# read fasta file
known <- readDNAStringSet(fasta)
# extract names
known.names <- names(known)
# strip names
kn.split <- str_split_fixed(known.names, ' ', 2)
known.names <- kn.split[,1]
head(known.names)

# read blastn file
blast.res <- read.blast(blastout)
head(blast.res)
blast.names <- unique(blast.res$query.id)
head(blast.names)

cat("Done reading files", '\n')

##------------Evaluate assembly against known starting sequences-----------

# get names of known genes that were mapped in assembly
captured <- which(known.names %in% blast.names)
length(captured)
captured.names <- known.names[captured]

# Compare length of assembled genes versus those missing
known.length <- width(known)
known.length.captured <- known.length[captured]
mean(known.length.captured)

known.length.missing <- known.length[-captured]
mean(known.length.missing)

# combine two dataframes for histogram
kcap <- data.frame(status="captured", length=known.length.captured)
kmis <- data.frame(status="missing", length=known.length.missing)

kdf <- rbind(kcap, kmis)


##-------------Report results-------------------------------------------

# plot
png("hist-missing-vs-captured-length.png")
ggplot(kdf, aes(length, fill=status)) +
    geom_density(alpha=0.2)
dev.off()

## Report results to file
cat('Genes', '\t', 'Count', '\t', 'Average length (base pairs)', '\n',
    'Starting', '\t', length(known.names), '\t', mean(width(known)), '\n',
    'Captured', '\t', length(captured), '\t', mean(known.length.captured), '\n',
    'Missing', '\t', length(known.length.missing), '\t', mean(known.length.missing), '\n',
    file="sim-assembly-results.txt")
