### CRISPR-Cas9 Finder Program

SYNOPSIS: A Perl program to find sgRNA sequences from Fasta or Bed format files.

  
##### sgRNA (single guide RNA)
The sgRNA is a 20 nt sequence used in CRISPR-Cas9 system used for targetting. 
a PAM sequence (NGG for Cas9) should be followed immedaitely after the sgRNA.

  
##### Strategy 
The searching for N20-NGG sequence is simple and straightforward.  
However evaluating the potential off-target effcts and the on-targeting activity is yet chanllenging. 
