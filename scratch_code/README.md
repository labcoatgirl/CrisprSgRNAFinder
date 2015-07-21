##### To Do List


1. Rewrite subroutines as individual PM files 
  * ~~read_fasta~~
  * ~~read_bed_file~~
  * ~~fetch_DNA_from_UCSC~~
  * ~~make a script and test the functions from the PM file~~
  * write detailed annoteation for each PM file
  * ~~DNAStuff module (complementary; DNA->RNA; and .... )~~
  
2. ~~Write the function to scan the sequence and return a list of sgRNA (20 nt) and its PAM~~
  * ~~Make the option to scan "plus" or "minus" or "both" strands~~
  * ~~Make the option to select top N sequences~~
  * ~~Calculate the GC contents (use GetGCPercentage function from DNAStuff module)~~

3. Write the function to evaluate the sgRNA sequence 
  * ~~How to update via the reference of the Array of Hash?~~
  * ~~write the module to Query GGGenome (without miss-match)~~ 
  * ~~calcualte the mismatch position by comparing with the original sequence~~
  * ~~tell the mismatch type (e.g: deletion, insertion, mismatch)~~
  * ~~return high-risk mismatch sequence number~~
  * ~~Create rigorous test data; (100bp-200bp will be fine); compare with other software~~
  * examine the high-risk mismatch sequence's location (e.g: exon; intron; intergenic)
  
