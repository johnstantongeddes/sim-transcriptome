sim-transcriptome
=================

Simulation of transcriptome assembly for gene expression analysis

The main script is `transcriptome-simulation.sh` that calls
`sample-fasta.R` to sample 100 sequences from the file ena.fasta. 
The script then simulates reads and runs a *standard* velvet-oases
transcriptome assembly pipeline, both with and without digital normalization.

Output files are

  sim-assembly-results-{suffix for each assembly}.txt 

which gives the number of known transcripts BLASTed by assembled transcripts and

  sim-assembly-results2-{suffix for each assembly}.txt 

which gives the mean length of known transcripts mapped, and the number of base pairs
in the assembled transcripts to the number bp in the known transcripts (roughly the 
number of isoforms with 1=one isoform)


The script `expression-simulation.sh` continues the above through mapping reads to 
the known and assembled transcripts to evaluate gene expression versus the simulated
expression levels from `rlsim`.

Programs required:
[rlsim](https://github.com/sbotond/rlsim)
[simNGS](http://www.ebi.ac.uk/goldman-srv/simNGS/)
[trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) and dependencies
[khmer](https://github.com/ctb/khmer)
[screed](https://github.com/ctb/screed)
[velvet](http://www.ebi.ac.uk/~zerbino/velvet/)
[oases](http://www.ebi.ac.uk/~zerbino/oases/)
[R](http://www.r-project.org/) and some Bioconductor packages
[cdhit](http://weizhong-lab.ucsd.edu/cd-hit/)

