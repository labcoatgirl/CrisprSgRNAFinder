#!/usr/bin/perl
# A perl program to search for CRISPR-Cas9 PAM sequence 
# could retrive sequences from UCSC DAS server
# supports FASTA file input
# query the off-targes number at super-fast gggenome engine (http://gggenome.dbcls.jp/)

use strict;
use warnings;
use 5.010;
use Getopt::Long;
#use LWP::UserAgent;
#use LWP::Simple;
#use XML::Simple;
#use JSON;

# use Pod::Usage;
use List::Util qw(min max);
use Term::ANSIColor;

our ($fasta_file, $bed_file, $genome, $mode,$length,$max_number, $help);

sub help 
{print <<USAGE;

--- Usage ---
	
perl cpfind.pl -b [bed file]| -f [fast file] -g [genome-build] -m [search mode] -l [sgRNA length] -n [top sequences]

--- Parameters ---

-b [the bed file] e.g: file1.bed
-f [the fasta file] e.g: file1.fa
-g [genome build] e.g: hg19|mm9|rn5 
-m [search mode] e.g: s|p| s: single mode| p: paired mode
-l [sgRNA length] e.g: 20 | default 20
-n [top sequences] e.g: 5 | return top 5 sgRNA sequences
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







# Expect a refernce of array of hash which returned by detectsgRNA
sub evaluatePAM{
	my $potential_target = shift;
	my $target_number = scalar (@$potential_target);
	for (my $i = 0; $i < $target_number ;$i++) {
		my $Num = $i + 1;
		
		my $perfect_match = GetCount($$potential_target[$i]->{sgRNA},0)-1 ;
		my $single_match  = GetCount($$potential_target[$i]->{sgRNA},1)-1 ;
		my $double_match  = GetCount($$potential_target[$i]->{sgRNA},2)-1 ;
		
		my $perfect_match_all = GetCount($$potential_target[$i]->{sgRNAandPAM},0)-1 ;
		my $single_match_all  = GetCount($$potential_target[$i]->{sgRNAandPAM},1)-1 ;
		my $double_match_all  = GetCount($$potential_target[$i]->{sgRNAandPAM},2)-1 ;
		my $triple_match_all  = GetCount($$potential_target[$i]->{sgRNAandPAM},3)-1 ;
		
		my $out_msg;
	    $out_msg .= "$Num of $target_number targets: ";
		$out_msg .= $$potential_target[$i]->{"sgRNAandPAM"};
		$out_msg .= " $$potential_target[$i]->{chr}:$$potential_target[$i]->{sta}-$$potential_target[$i]->{end}\n";
		$out_msg .= "sgRNA without PAM missmatch number (0-1-2): $perfect_match $single_match $double_match \n";
		$out_msg .= "sgRNA with  PAM missmatch number (0-1-2-3): $perfect_match_all $single_match_all $double_match_all $triple_match_all\n";
		print $out_msg;
	}
}


