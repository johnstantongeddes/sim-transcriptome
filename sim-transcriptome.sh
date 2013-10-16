#!/bin/bash

##################################################################################################
## Script to generate simulate Illumina RNAseq data and validate transcriptome assembly
## Author: John Stanton-Geddes
## Created: 2013-08-12
###################################################################################################

# Set PYTHONPATH to khmer module
 
export PYTHONPATH=/opt/software/khmer/python

##--------------------Simulate data------------------------------------------------

# Fasta file of known mRNA transcripts for simulation
# In this example, I use a file with 1000 *Arabidopsis* mRNA transcripts downloaded from [European Nucleotide Archive](http://www.ebi.ac.uk/ena/home)

knownfasta="ena.fasta"

# Randomly select X sequences; output "known.fasta"
Rscript sample-fasta.R $knownfasta 100

## Trim fasta titles to just ID
cut -f1 -d" " known.fasta > known1.fasta

## Assign expression-level to each sequence using the `sel` tool in [rlsim](https://github.com/sbotond/rlsim)

sel -d "1.0:g:(1500, 3)" known1.fasta > known-sim.fasta
rm known.fasta known1.fasta # Clean-up temp files
echo -e "Done with expression levels: " `date` "\n"

## Simulate fragments created during library prep using `rlsim`. Based on length of fragments estimated when using FLASH to pair reads, use empirical length distribution of 180 bp (SD: 20bp).

rlsim -v -n 1500000 -d "1:n:(180, 20, 100, 500)" known-sim.fasta > known-sim-frags.fasta
echo -e "Done with simulated fragmentation levels: " `date` "\n"



##-----------------Assemble simulated data---------------------------------------------

## Generate simulated Illumina paired-end reads using [simNGS](http://www.ebi.ac.uk/goldman-srv/simNGS/) and a default runfile provided with simNGS
# output is 'reads_end1.fq' and 'reads_end2.fq'
cat known-sim-frags.fasta | simNGS -p paired -o fastq -O reads /opt/software/simNGS/data/s_3_4x.runfile
echo -e "Done with Illumina read simulation: " `date` "\n"

## Run through TrimGalore for adapter cutting and quality trimming
trim_galore --quality 20 --phred33 --fastqc --length 20 --paired reads_end1.fq reads_end2.fq
echo -e "Done with quality trimming: " `date` "\n"

## Run digital normalization with khmer

# interleave file
python /opt/software/khmer/scripts/interleave-reads.py reads_end1_val_1.fq reads_end2_val_2.fq > sim-interleaved.fastq

# run diginorm specifying paired-end reads (-p) with coverage threshold and kmer of 20
python /opt/software/khmer/scripts/normalize-by-median.py -R diginorm-final.out -p -C 20 -k 20 -N 4 -x 4e9 sim-interleaved.fastq
echo -e "Done with interleaving files:" `date` "\n"

##----------Assemble transcriptome with velvet-oases------------------------

# Assemble for all reads
# velveth
velveth sim-oases-21 21 -fastq -interleaved -shortPaired sim-interleaved.fastq
# velvetg
velvetg sim-oases-21 -read_trkg yes -exp_cov auto
# oases
oases sim-oases-21 -ins_length 180
# check basic stats
python /opt/software/khmer/sandbox/assemstats2.py 100 sim-oases-21/transcripts.fa > oases-21-assem-stats.txt

echo -e "Done with velvet-oases for normalized reads:" `date` "\n"

# Assemble for digitally-normalized reads
# velveth
velveth sim-oases-norm-21 21 -fastq -interleaved -shortPaired sim-interleaved.fastq.keep
# velvetg
velvetg sim-oases-norm-21 -read_trkg yes -exp_cov auto
# oases
oases sim-oases-norm-21 -ins_length 180
# check basic stats
python /opt/software/khmer/sandbox/assemstats2.py 100 sim-oases-norm-21/transcripts.fa > oases-norm-21-assem-stats.txt

echo -e "Done with velvet-oases for all reads:" `date` "\n"


##----------------BLAST known starting sequences against assembly---------------

## BLAST assembled transcripts against true transcripts
# Make blast database
makeblastdb -dbtype nucl -in sim-oases-21/transcripts.fa
# specify output in tabular format using `-outfmt 6` so that file can be read by RFLPtools in R
# for post-processing
blastn -query known-sim.fasta -db sim-oases-21/transcripts.fa -outfmt 6 -out blast-oases-21.txt
blastn -query known-sim.fasta -db sim-oases-21/transcripts.fa -out blast-oases-21.html

## Repeat for normalized assembly
makeblastdb -dbtype nucl -in sim-oases-norm-21/transcripts.fa
blastn -query known-sim.fasta -db sim-oases-norm-21/transcripts.fa -outfmt 6 -out blast-oases-norm-21.txt
blastn -query known-sim.fasta -db sim-oases-norm-21/transcripts.fa -out blast-oases-norm-21.html

 
# R script to extract top assembled transcript against each starting sequence
# input: (1) known fasta file (2) blast results
# output: 
#   sim-assembly-results1-{blast file prefix}.txt  - number of starting transcripts that are recaptured and length of captured versus missing transcripts
#   sim-assembly-results2-{blast file prefix}.txt  - mean proportion of original transcript mapped and mean bp mapped to starting number bp
Rscript sim-assembly-eval.R known-sim.fasta blast-oases-21.txt
Rscript sim-assembly-eval.R known-sim.fasta blast-oases-norm-21.txt


##---------------Map reads against against known transcripts----------------------------

## Map simulated reads to known transcripts using TopHat-Cufflinks
# Create index with Bowtie
bowtie2-build known-sim.fasta bowtie-known
# Run tophat to map reads to reference
tophat -o tophat_known_out bowtie-known reads_end1_val_1.fq reads_end2_val_2.fq 
# Run cufflinks for transcript discovery
cufflinks -o cufflinks_known_out tophat_known_out/accepted_hits.bam


## Map simulated reads to velvet-oases assembled transcripts using TopHat-Cufflinks
# Create index with Bowtie
bowtie2-build sim-oases-21/transcripts.fa bowtie-assembled
# Run tophat to map reads to reference
tophat -o tophat_assembled_out bowtie-assembled reads_end1_val_1.fq reads_end2_val_2.fq 
# Run cufflinks for transcript discovery
cufflinks -o cufflinks_assembled_out tophat_assembled_out/accepted_hits.bam
## Hmmm...cufflinks only reports counts for ~30 of 100 genes/isoforms...where are the rest?

## Map simulated reads to velvet-oases assembled transcripts using TopHat-Cufflinks
# Create index with Bowtie
bowtie2-build sim-oases-norm-21/transcripts.fa bowtie-assembled-norm
# Run tophat to map reads to reference
tophat -o tophat_assembled_norm_out bowtie-assembled-norm reads_end1_val_1.fq reads_end2_val_2.fq 
# Run cufflinks for transcript discovery
cufflinks -o cufflinks_assembled_norm_out tophat_assembled_norm_out/accepted_hits.bam
## Hmmm...cufflinks only reports counts for ~30 of 100 genes/isoforms...where are the rest?

## Possibly due to multiple alignment of reads to assembled transcripts incorrectly inferring multiple isoforms...

cat tophat_known_out/align_summary.txt
# 6.5% aligned pairs have multiple alignments

cat tophat_assembled_out/align_summary.txt
# 52.7% aligned pairs have multiple alignments

cat tophat_assembled_norm_out/align_summary.txt
# 59.4% aligned pairs have multiple alignments







##---------------------- Map with BWA------------------------------------

# Build index with BWA
bwa index known-sim.fasta

# Map paired reads
bwa mem known-sim.fasta reads_end1_val_1.fq reads_end2_val_2.fq > BWA-known-mapped-pe.sam

# Convert to BAM
samtools faidx known-sim.fasta # index
samtools import known-sim.fai BWA-known-mapped-pe.sm sim-mapped-pe.bam # sam -> bam
samtools sort sim-mapped-pe.bam sim-mapped-pe.sorted # sort BAM
samtools index sim-mapped-pe.sorted.bam # index

# convert known fasta file into BED using script from http://ged.msu.edu/angus/tutorials-2013/rnaseq_bwa_counting.html?highlight=bwa
python /home/projects/climate-cascade/scripts/make_bed_from_fasta.py known-sim.fasta > knownfasta.bed

# count reads that have mapping quality of 30 or better `-q 30` to the known transcripts
multiBamCov -q 30 -p -bams sim-mapped-pe.sorted.bam -bed knownfasta.bed > BWA-counts-known.txt

# Calculate Transcripts per million (TPM) (Wagner et al. 2012 Theory. Biosci) 
Rscript TPM.R BWA-sim-counts.txt

##---------------Map reads against against assembled transcripts----------------------------

# Index
bwa index sim-oases-21/transcripts.fa

# Map paired reads
bwa mem sim-oases-21/transcripts.fa reads_end1_val_1.fq reads_end2_val_2.fq > oases-mapped-pe.sam

# Convert to BAM
samtools faidx sim-oases-21/transcripts.fa # index
samtools import sim-oases-21/transcripts.fai oases-mapped-pe.sam oases-mapped-pe.bam # sam -> bam
samtools sort oases-mapped-pe.bam oases-mapped-pe.sorted # sort BAM
samtools index oases-mapped-pe.sorted.bam # index

# convert oases transcripts to bed 
python /home/projects/climate-cascade/scripts/make_bed_from_fasta.py sim-oases-21/transcripts.fa > sim-oases-21-transcripts.bed

# count reads that have mapping quality of 30 or better to the assembled transcripts
multiBamCov -q 30 -p -bams oases-mapped-pe.sorted.bam -bed sim-oases-21-transcripts.bed > BWA-counts-oases.txt


# Calculate Transcripts per million (TPM) (Wagner et al. 2012 Theory. Biosci) 
Rscript TPM.R BWA-counts-oases.txt


## This results in reads being mapped to 203 transcripts as oases identifies multiple isoforms 
## for genes that do not actually exist
## Reduce oases assembly by selecting longest transcript for each locus 
Rscript oases-reduce.r sim-oases-21/transcripts.fa blast-oases-21.txt 
# output file: oases-transcripts-kept.fa


##----------Map reads against unique assembled transcripts with BWA----------------
# Index
bwa index oases-transcripts-kept.fa
# Map paired reads
bwa mem oases-transcripts-kept.fa reads_end1_val_1.fq reads_end2_val_2.fq > oases-transcripts-kept-mapped.sam
# Convert to BAM
samtools faidx oases-transcripts-kept.fa # index
samtools import oases-transcripts-kept.fai oases-transcripts-kept-mapped.sam oases-transcripts-kept-mapped.bam # sam -> bam
samtools sort oases-transcripts-kept-mapped.bam oases-transcripts-kept-mapped-sorted # sort BAM
samtools index oases-transcripts-kept-mapped-sorted.bam # index
# convert oases transcripts to bed 
python /home/projects/climate-cascade/scripts/make_bed_from_fasta.py oases-transcripts-kept.fa > oases-transcripts-kept.bed
# count reads that have mapping quality of 30 or better to the assembled transcripts
multiBamCov -q 30 -p -bams oases-transcripts-kept-mapped-sorted.bam -bed oases-transcripts-kept.bed > BWA-counts-kept.txt

## Compare read counts for mapping to known vs assembled transcripts
# first extract counts for the 78 known loci that correspond to the 78 that were assembled 
##### Rscripts extract-known.R 

##### Rscript known-vs-assembled.R BWA-sim-counts.txt BWA-counts-kept.txt
