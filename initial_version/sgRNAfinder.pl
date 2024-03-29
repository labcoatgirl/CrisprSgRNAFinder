#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
use Getopt::Long;
# use Pod::Usage;
use List::Util qw(min max);
use Term::ANSIColor;


our ($fasta_file, $bed_file, $genome, $mode,$length,$max_number, $help, $output);

sub help 
{print <<USAGE;

Usage:

perl sgRNAfinder.pl -b [bed file]| -f [fasta file] -o [output file] -g [genome-build] -m [search mode] -l [sgRNA length] -n [top sequences] 

Parameters:

-b [the input bed file] 				e.g: example.bed
-f [the input fasta file]				e.g: example.fa
-o [the output file]					e.g: example_out.txt
-g [genome build] 					[default: hg19] | e.g: hg19;mm9;rn5; 
-m [search mode]					[default: s]	| e.g: [s|p] s: single mode; p: paired mode	
-l [sgRNA guide sequence length				[default: 20]	| the default length of guide sequence is 20 nt  
-n [top sequences] 					[default: 10]	| return top 10 sgRNA sequences 

-h [help]

USAGE
exit;
}

GetOptions (			
"fasta|f=s"		=> \$fasta_file,
"bed|b=s"		=> \$bed_file,
"genome|g=s"	=> \$genome,
"mode|m=s"		=> \$mode,
"length|l=i"	=> \$length,
"max|n=s"		=> \$max_number,
"help|h"		=> \&help,
"output|o=s"	=> \$output,
) or die ("\nPlease check the arguments; perl cpfind.pl -h for help\n\n");




defined ($fasta_file or $bed_file or $genome or $mode or $length or $max_number or $help) or print &help;
if (!$fasta_file and !$bed_file){
	print "Error: please provide the input file (Bed or Fasta format)\n\n"
} 



# ==================================
# = Initialize default Parameters  =
# ==================================
$genome     = defined($genome) 		?	$genome 	: 	"hg19";
$mode       = defined($mode) 		?	$mode 		: 	"s";
$length     = defined($length) 		?	$length		:	20;
$max_number = defined($max_number) 	?	$max_number :	10;


# ==================================
# =   Read Fasta File if provided  =
# ==================================

if (defined $fasta_file){
	print colored("Now is reading $fasta_file\n", "bold blue");
	print colored("Sequence is based on $genome\n", "bold blue");
	my $seq = readfasta($fasta_file);
	my $seq_num = scalar(@$seq);
	print colored("Finish reading $fasta_file\n$seq_num sequence(s) have been read\n", "bold blue");
	
	for my $item (@$seq) 
	{
		my $mytargets = detectsgRNA($item);
		evaluatePAM($mytargets);
		print "=========================================\n";
	}

}


# ==================================
# =    Read Bed File if provided   =
# ==================================

if (defined $bed_file){
	print colored("Now is reading $bed_file\n", "bold blue");
	print colored("Sequence is based on $genome\n", "bold blue");
	my $seq = readbed($bed_file);
	my $seq_num = scalar(@$seq);
	print colored("Finish reading $bed_file\n$seq_num sequence have been read\n", "bold blue");
	for my $item (@$seq) 
	{	
		my $mytargets = detectsgRNA($item);
		evaluatePAM($mytargets);
		print "=========================================\n";
	}
	
}





