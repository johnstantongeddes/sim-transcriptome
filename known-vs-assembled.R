## Script to compare counts of reads mapped to known versus oases assembled transcripts

##--------Read files----------------

known.counts <- read.table("sim-transcriptome-counts-BWA.txt", header=FALSE, sep="\t")
colnames(known.counts) <- c("transcript", "start", "stop", "count")
head(known.counts)
dim(known.counts)
str(known.counts)

assembled.counts <- read.table("BWA-counts-oases.txt", header=FALSE, sep="\t")
colnames(assembled.counts) <- c("transcript", "start", "stop", "count")
head(assembled.counts)
dim(assembled.counts)

