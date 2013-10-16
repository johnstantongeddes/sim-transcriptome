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
dim(blast.res)
blast.names <- unique(blast.res$query.id)
head(blast.names)
length(unique(blast.names))

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


##--------Evaluate quality of assembly-------------------------------
# of known transcripts that are mapped, what proportion of their length is mapped?
# does Oases infer isoforms that don't really exist?

# Dataframe to collect results
qual.df <- matrix(ncol=3, nrow=0)
colnames(qual.df) <- c("gene", "length.mapped", "bp.mapped")

for(i in 1:length(unique(blast.res$query.id))) {
    # get number of transcripts
    bi <- blast.res$query.id[i]
    bin <- blast.res[blast.res$query.id==bi,]
    nrow(bin)
    # length of known transcript
    kl <- width(known[which(known.names==bi)])
    # what proportion of known transcript mapped by assembled transcripts
    prop.mapped <- round((max(bin$q.end) - min(bin$q.start))/kl,2)
    # how many base pairs of assembled transcript relative to known (e.g. incorrect isoforms)?
    bp.mapped <- sum(bin$alignment.length)/kl
    # report
    qual.df <- rbind(qual.df, c(bi, prop.mapped, bp.mapped))
}

qual.df <- as.data.frame(qual.df)
qual.df$length.mapped <- as.numeric(as.character(qual.df$length.mapped))
qual.df$bp.mapped <- as.numeric(as.character(qual.df$bp.mapped))
str(qual.df)

# Mean length of original transcript mapped
mean(qual.df$length.mapped)
# Mean copies of original transcript mapped. 1=mapped once, <1 mapped to less than full length,
#  > 1 mapped multiple times (e.g 2=mapped twice)
mean(qual.df$bp.mapped)

##-------------Report results-------------------------------------------

# prefix for output files
pre <- str_split_fixed(blastout, "\\.", 2)[,1]

# plot
png(paste("hist-", pre, "-missing-vs-captured-length.png", sep=""))
ggplot(kdf, aes(length, fill=status)) +
    geom_density(alpha=0.2)
dev.off()

## Report results to file
cat('Transcripts', '\t', 'Count', '\t', 'Average length (base pairs)', '\n',
    'Starting', '\t', length(known.names), '\t', mean(width(known)), '\n',
    'Captured', '\t', length(captured), '\t', mean(known.length.captured), '\n',
    'Missing', '\t', length(known.length.missing), '\t', mean(known.length.missing), '\n',
    file=paste("sim-assembly-results-", pre, ".txt", sep=""))

cat('Mean proportion of original transcript mapped', '\t', 'Proportion bp assembled to starting bp', '\n',
    round(mean(qual.df$length.mapped),2), '\t', round(mean(qual.df$bp.mapped),2), '\n',
    file=paste("sim-assembly-results2-", pre, ".txt", sep=""))


