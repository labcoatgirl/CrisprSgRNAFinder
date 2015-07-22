#!/usr/bin/perl
# A perl program to search for CRISPR-Cas9 PAM sequence 
# could retrive sequences from UCSC DAS server
# supports FASTA file input
# query the off-targes number at super-fast gggenome engine (http://gggenome.dbcls.jp/)

use strict;
use warnings;
use 5.010;
use Getopt::Long;
# use Pod::Usage;
use List::Util qw(min max);
use Term::ANSIColor;

our ($fasta_file, $bed_file, $genome, $mode,$length,$max_number, $help);

sub help 
{print <<USAGE;

--- Usage ---
perl sgRNAfinder.pl -b [bed file]| -f [fasta file] -g [genome-build] -m [search mode] -l [sgRNA length] -n [top sequences]

--- Parameters ---
-b [the bed file] 					e.g: example.bed
-f [the fasta file]					e.g: example.fa
-g [genome build] 					e.g: hg19;mm9;rn5;  
-m [search mode] 					e.g: s|p| s: single mode| p: paired mode
-l [sgRNA guide sequence length]			e.g: 20 | default 20
-n [top sequences] 					e.g: 5 	| return top 5 sgRNA sequences
-o [output file]					e.g: example_out.txt
-h [help]

USAGE
exit;
}

GetOptions (			
"fasta|f=s" => \$fasta_file,
"bed|b=s" => \$bed_file,
"genome|g=s" => \$genome,
"mode|m=s" => \$mode,
"length|l=i" => \$length,
"max|n=s" => \$max_number,
"help|h" => \&help,
) or die ("\nPlease check the arguments; perl cpfind.pl -h for help\n\n");

defined ($fasta_file or $bed_file or $genome or $mode or $length or $max_number or $help) or &help;
if (!$fasta_file and !$bed_file){print "\nPlease provide the input file (bed or Fasta format)\n\n"} 



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





