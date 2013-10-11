## Script to pull longest transcript for oases for each locus when comparing against known loci used for simulation


##-------------Load libraries--------------------------------
library(Biostrings)

##-------------Get command-line arguments and read files------------
args <- commandArgs(trailingOnly = TRUE)
fasta <- args[1]
blastout <- args[2]
transcripts <- args[3]


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


head(blast.res)

blast.sub <- blast.res[0,]

for(i in unique(blast.res$query.id)) {
    # get all transcripts that match to single query
    bin <- blast.res[blast.res$query.id==i,]
    # extract longest transcript, or simply first if multiple of equal length
    longest <- bin[which(bin$alignment.length == max(bin$alignment.length))[1], ]
    blast.sub <- rbind(blast.sub, longest)
}

if(nrow(blast.sub) != length(unique(blast.res$query.id))) stop("missing transcripts")

# Read Oases transcripts

# read fasta file of oases transcripts
oases <- readDNAStringSet(transcripts)

# extract names
oases.names <- names(oases)


# Pull out oases transcripts from longest BLAST matches and save to file
for(j in 1:nrow(blast.sub)) {
    oases.select <- oases[which(names(oases) == blast.sub$subject.id[j]) ,]
    writeXStringSet(oases.select, file="oases-transcripts-kept.fa", append=TRUE)
}



    


     
