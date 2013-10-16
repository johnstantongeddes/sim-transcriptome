## Script to pull longest transcript for oases for each locus when comparing against known loci used for simulation


##-------------Load libraries--------------------------------
library(Biostrings)
library(stringr)
library(RFLPtools)

##-------------Get command-line arguments and read files------------
args <- commandArgs(trailingOnly = TRUE)
fasta <- args[1]
blastout <- args[2]
transcripts <- args[3]

## Read known fasta file
known <- readDNAStringSet(fasta)
# extract names
known.names <- names(known)
# strip names
kn.split <- str_split_fixed(known.names, ' ', 2)
known.names <- kn.split[,1]
head(known.names)

## Read Oases transcripts
# read fasta file of oases transcripts
oases <- readDNAStringSet(transcripts)
head(oases)
# extract names
oases.names <- names(oases)

## Read blastn file of oases transcripts mapped to known transcripts
blast.res <- read.blast(blastout)
head(blast.res)
blast.names <- unique(blast.res$query.id)
head(blast.names)



cat("Done reading files", '\n')

##------Extract longest transcript from each locus--------------------
blast.sub <- blast.res[0,]

for(i in unique(blast.res$query.id)) {
    # get all transcripts that match to single query
    bin <- blast.res[blast.res$query.id==i,]
    # extract longest transcript, or simply first if multiple of equal length
    longest <- bin[which(bin$alignment.length == max(bin$alignment.length))[1], ]
    blast.sub <- rbind(blast.sub, longest)
}

if(nrow(blast.sub) != length(unique(blast.res$query.id))) stop("missing transcripts")

dim(blast.sub)
# By BLAST, 86 of original 100 transcripts are recaptured in assembled transcripts 
length(unique(blast.sub$subject.id))
# but these are represented by only 78 unique assembled loci 
# visually checked and duplicates are closely related known transcripts
# e.g. 


# Extract out unique oases transcripts from longest BLAST matches and save to file
for(j in unique(blast.sub$subject.id)) {
    oases.select <- oases[which(names(oases) == j) ,]
    writeXStringSet(oases.select, file="oases-transcripts-kept.fa", append=TRUE)
}


# Extract names of known transcripts that correspond to the 'oases-transcripts-kept.fa'

# read file
kept <- read.table("oases-transcripts-kept.fa", header=FALSE, sep="\t")
head(kept)

head(oases.select)



    


     
