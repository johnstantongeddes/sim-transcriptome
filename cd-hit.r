## Pull longest transcript for oases for each locus

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
oases <- readDNAStringSet("sim-oases-21/transcripts.fa")

# extract names
oases.names <- names(oases)


# Pull out oases transcripts from longest BLAST matches and save to file
for(j in 1:nrow(blast.sub)) {
    oases.select <- oases[which(names(oases) == blast.sub$subject.id[j]) ,]
    writeXStringSet(oases.select, file="oases-transcripts-kept.fa", append=TRUE)
}
# Export


    


     
