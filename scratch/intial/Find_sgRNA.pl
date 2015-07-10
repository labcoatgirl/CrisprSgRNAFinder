#!/usr/bin/perl
# A perl program to search for CRISPR-Cas9 PAM sequence 
# could retrive sequences from UCSC DAS server
# supports FASTA file input
# query the off-targes number at super-fast gggenome engine (http://gggenome.dbcls.jp/)

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use LWP::UserAgent;
use LWP::Simple;
use XML::Simple;
use JSON;

# use Pod::Usage;
use List::Util qw( min max );
use Term::ANSIColor;

our ($fasta_file, $bed_file, $genome, $mode,$length,$max_number, $help);

sub help 
{
print <<USAGE;

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
if (!$fasta_file and !$bed_file){print "\nPlease provide the input file - bed file or fasta file\n\n"} 



# ==================================
# = Initialize default Parameters  =
# ==================================
$genome     = defined($genome) ? $genome : "hg19";
$mode       = defined($mode) ? $mode : "s";
$length     = defined($length) ? $length:20;
$max_number = defined($max_number) ? $max_number :10;


# ==================================
# =   Read Fasta File if provided  =
# ==================================

if (defined $fasta_file){
	print colored("Now is reading $fasta_file\n", "bold blue");
	print colored("Sequence is based on $genome\n", "bold blue");
	my $seq = readfasta($fasta_file);
	my $seq_num = scalar(@$seq);
	print colored("Finish reading $fasta_file\n$seq_num sequence have been read\n", "bold blue");
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


sub readfasta 
{
    (my $file)= shift;
	my %sequence;
	my $header="";
	my $temp_seq;
	
	#suppose fasta files contains multiple sequences;
	 
	open (IN, "<$file") or die "couldn't open the file $file $!";
	while (<IN>)
	{	
		$_ =~ s/\R//g;
		next if /^\s*$/; #skip empty line 
		if ($_ =~ /^>/)  #when see head line
		{	
			$header= $_;
			if ($sequence{$header}){print colored("#CAUTION: SAME FASTA HAS BEEN READ MULTIPLE TIMES.\n#CAUTION: PLEASE CHECK FASTA SEQUENCE:$header\n","red")};
			if ($temp_seq) {$temp_seq=""} # If there is alreay sequence in temp_seq, empty the sequence file
			
		}
		else # when see the sequence line 
		{
		   s/\s+//g;
		   $temp_seq .= $_;
		   $sequence{$header}=$temp_seq; #update the contents
		}
	
	}
	
	# Reformat the sequences and make an array of array format 
	my @AoA;
	foreach my $k (sort keys %sequence) {
		my ($chr,$sta,$end);
		if ($k =~ /.*(chr\w+):(\d+)-(\d+)/g ){
					$chr = $1 ;
					$sta = $2 ;
					$end = $3 ;
			}else{
					$chr = "0" ;
					$sta = "0" ;
					$end = "0" ;
			}
#		print "$sequence{$k} $chr $sta $end \n"; 
		push @AoA, [$sequence{$k},$chr,$sta,$end] 
		
	}
	
	undef %sequence;
	
	return \@AoA;
}


sub Parse_UCSC
{
	my ($chr,$sta,$end) = @_;
	my ($link) ="http://genome.ucsc.edu/cgi-bin/das/"."$genome"."/dna?segment=$chr:$sta,$end";
	my $xml = XML::Simple->new();
	my $data = $xml->XMLin(get($link));
	my $sequence = $data-> {SEQUENCE}->{DNA}->{content};
	$sequence =~ s/(\s+)//gi;
	$sequence = uc($sequence);
	return $sequence;
}

sub readbed 
{
	my $file = shift;
	my @AoA;
	open (IN, "<$file") or die "couldn't open the file $file $!";
	while (<IN>)
	{	
		#A typical bed format file looks like this 
		#track name=pairedReads description="Clone Paired Reads" useScore=1
		#chr22 1000000 1000500 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
		#chr22 1000600 1000800 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
		
		chop;
		next if /^\s*$/;  #skip empty line 
		next unless /^chr.*/; 
		my @bedcolumns = split(/\s+/, $_);
#		print $bedcolumns[0],$bedcolumns[1],$bedcolumns[2];
		my $seq= Parse_UCSC($bedcolumns[0],$bedcolumns[1],$bedcolumns[2]);
		push @AoA, [$seq,$bedcolumns[0],$bedcolumns[1],$bedcolumns[2]] 
#		print $bedcolumns[0];
	}
	
	return \@AoA;
}


# ================================
# ==== Detect PAM sequences =====
# ================================
# detectPAM receive sequence and chromosome positions 
# and returns an array of hash

sub complementary {
  my ($sequence) = shift;
  my $reverse_sequence = reverse($sequence); 
  $reverse_sequence =~ tr/ACGTRYNX/TGCAYRNX/;
  return $reverse_sequence;
}


sub detectsgRNA
{
	my $AoA = shift;
	my ($seq,$chr,$chr_sta,$chr_end) = @$AoA;
	my @targets;
	
	unless (defined $chr_sta)
	{
	$chr_sta = 0; # use 0 based system 
	$chr = 0;		
	}
	
	while ($seq =~ /(?=gg)/gi) #First process + strand 
	{  	
		# An illustration for sgRNA position 
		# positions and sequences are as shown below:
		# if start position is not provided, use 0 based system
		# 0 - 1 - 2 - 3 - 4 - 5 - 6 - 7 
		# A - C - C - T - G - G - C - C
		# pos($seq) = 4 -> the G-G starts at position 4
		# suppose sgRNA lenght is 2 bp not including NGG (usually we have 20 nt for sequence before GG) 
		# crispr_sta should start at position 1
		# crispr_sta = 4 - length(2) - 1   
		
		# A general form is  
		# crispr_sta = 0 + pos(seq) - length(sgRNA) - 1 
		# crispr_end = 0 + pos(seq) + 1 
		
		
		# In a different situation where the start position is not 0
		# 1306 1307 1308 1309 1310 1311 1312 1313 1314 1315
		#  A    C     C    T   G     C    G    G    A    A
		#  0    1     2    3   4     5    6    7    8    9
		# pos($seq) = 6 -> the G-G starts at position 6
		# suppose sgRNA lenght is 2 bp not including NGG 
		# crispr_sta should start at position 1309
		# crispr_sta = 6 - length(2) - 1 + start = 6 - 2 - 1 + 1306 = 1309 
		# crispr_sta = start + pos(seq) - length(sgRNA) - 1 
		# crispr_end = start + pos(seq) + 1 
		
		my $crispr_sta = $chr_sta + pos($seq) - $length - 1;
		my $crispr_end = $chr_sta + pos($seq) + 1;
		
		# double check the total length
		# total_length = $crispr_end - $crispr_sta + 1 =  3 + length (OK!)
		
		
       	# when estimate PAM sequences, we do not use whole sequences (sgRNA + PAM)
		# but rather check the sequence before NGG (sgRNA)
		# The total sequnce is saved as sgRNAandPAM
		
		# Continue our illustration above
		# The sgRNA should be the relative position to the start position
		# sg_start = pos($seq) - $length - 1
		# sg_end   = pos($seq) - $length - 1 + length - 1 = pos($seq) - 2
		
		my $sg_start =  pos($seq) - $length - 1; 
		my $sg_end   =  pos($seq) - 2;
		
		if ($crispr_sta >= $chr_sta)
		{

			my $target=substr($seq,$sg_start,$length);
			my $full_target =substr($seq,$sg_start,$length+3);
			my $temp_target;			
			$temp_target ->{"chr"} = $chr;
			$temp_target ->{"sta"} = $crispr_sta;
			$temp_target ->{"end"} = $crispr_end;
			$temp_target ->{"sgRNA"} = $target;
			$temp_target ->{"sgRNAandPAM"} = $full_target;
			$temp_target ->{"strand"} = "+";
			push @targets, $temp_target;	# Create Array of Hash
		} 		
	}	
	
	while ($seq =~ /(?=cc)/gi) #Next process the minus strand
	{  	
		# 12 - 13 - 14 -15-16- 17 -18 - 19
		# 0  - 1  - 2 - 3 - 4 - 5 - 6 - 7 
		# A  - C  - C - T - G - G - C - C
		# First CC start at position 1 - which is also the end position of GG
		
		my $crispr_sta = $chr_sta + pos($seq);
		my $crispr_end = $chr_sta + pos($seq)+ 2 + $length;
		
		# suppose the sgRNA length is 2
		# crispr_end = 12 + 1 + 2 + 2 = 17 (OK)
		
		my $fetch_start =  pos($seq)+3;
		my $fetch_end   =  pos($seq)+ $length + 2;
		
		if ($fetch_end <= length($seq) -1 )
		{
			my $target=complementary(substr($seq,$fetch_start,$length)); 			# Get the complementary sequence
			my $full_target =complementary(substr($seq,pos($seq),$length+3));       # Get the complementary sequence
			my $temp_target ={};
			$temp_target -> {"chr"} = $chr;
			$temp_target -> {"sta"} = $crispr_sta;
			$temp_target -> {"end"} = $crispr_end;
			$temp_target -> {"sgRNA"} = $target;
			$temp_target -> {"sgRNAandPAM"} = $full_target;
			$temp_target -> {"strand"} = "-";
			push @targets, $temp_target;	
		} 		
	}	
	return \@targets;
}


sub GetCount
{
	my ($target_sequence,$mismatch) = @_;
	my $url_base = 'http://gggenome.dbcls.jp/'.$genome;
	my $perfect_url;
	if ($mismatch eq 0){
		$perfect_url=$url_base.'/'.$target_sequence.'.json'
	}else{
		$perfect_url=$url_base.'/'.$mismatch.'/'.$target_sequence.'.json';
	}
	
	my $perfect_ua = LWP::UserAgent->new();
	my $perfect_response = $perfect_ua->get($perfect_url);
	my $perfect_json = JSON->new;
	my $perfect_data = $perfect_json->decode($perfect_response->content());
	my $total_count=$perfect_data->{'summary'}[0]{'count'}+$perfect_data->{'summary'}[1]{'count'};
	return $total_count;
	
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


