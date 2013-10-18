#!/bin/bash

##################################################################################################
## Script to generate simulate Illumina RNAseq data and validate transcriptome assembly
## Author: John Stanton-Geddes
## Created: 2013-08-12
###################################################################################################

# Set PYTHONPATH to khmer module
 
export PYTHONPATH=/opt/software/khmer/python

# Make directory for intermediate results
mkdir -p out

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

## Generate simulated Illumina paired-end reads using [simNGS](http://www.ebi.ac.uk/goldman-srv/simNGS/) and a default runfile provided with simNGS
# output is 'reads_end1.fq' and 'reads_end2.fq'
cat known-sim-frags.fasta | simNGS -p paired -o fastq -O reads /opt/software/simNGS/data/s_3_4x.runfile
echo -e "Done with Illumina read simulation: " `date` "\n"

# move intermediate files
mv rlsim_report.json out
mv sel_report.pdf out

##---------------QC and diginorm simulated Illumina reads----------------------

## Run through TrimGalore for adapter cutting and quality trimming
trim_galore --quality 20 --phred33 --fastqc --length 20 --paired reads_end1.fq reads_end2.fq
echo -e "Done with quality trimming: " `date` "\n"

# move intermediate files
mv reads_end*.fq out
mv reads_end*.fq_trimming_report.txt out
mv reads_end*_val_*.fq_fastqc out
mv reads_end*_val_*.fq_fastqc.zip out


## Run digital normalization with khmer

# interleave file
python /opt/software/khmer/scripts/interleave-reads.py reads_end1_val_1.fq reads_end2_val_2.fq > sim-interleaved.fastq

# run diginorm specifying paired-end reads (-p) with coverage threshold and kmer of 20
python /opt/software/khmer/scripts/normalize-by-median.py -R diginorm-final.out -p -C 20 -k 20 -N 4 -x 4e9 sim-interleaved.fastq
echo -e "Done with interleaving files:" `date` "\n"

mv diginorm-final.out out 


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


##-----------CD-HIT on assembled transcripts--------------------------------------------

cd-hit-est -i sim-oases-21/transcripts.fa -o sim-oases-21-cdhit.fa -c 0.95 -n 8
cd-hit-est -i sim-oases-norm-21/transcripts.fa -o sim-oases-norm-21-cdhit.fa -c 0.95 -n 8

## BLAST reduced transcripts to known assembly
makeblastdb -dbtype nucl -in sim-oases-21-cdhit.fa
blastn -query known-sim.fasta -db sim-oases-21-cdhit.fa -outfmt 6 -out blast-oases-21-cdhit.txt
blastn -query known-sim.fasta -db sim-oases-21-cdhit.fa -out blast-oases-21-cdhit.html

Rscript sim-assembly-eval.R known-sim.fasta blast-oases-21-cdhit.txt
