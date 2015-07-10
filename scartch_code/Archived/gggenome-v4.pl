#!/usr/bin/perl
# A local application to design the Crispr sgRNA sequence 
# the off-targets sites were serached via gggenome engine via REST API
# the candidate targets were scored according to Feng Zhang lab's formula

use Getopt::Long;
use strict;
use warnings;
use LWP::UserAgent;
use JSON;
use List::Util qw( min max );

my $length="";
my $fasta_file="";
my $bed_file="";
my $species="hg19";

GetOptions ("length|l=i" => \$length,
			"species|s=s" => \$species,
			"fasta|f=s" => \$fasta_file,
			"bed|b=s" => \$bed_file,
			);

my %table =
(
	'1' => 0,
	'2' => 0,
	'3' => 0.014,
	'4' => 0,
	'5' => 0,
	'6' => 0.395,
	'7' => 0.317,
	'8' => 0,
	'9' => 0.389,
	'10' => 0.079,
	'11' => 0.445,
	'12' => 0.508,
	'13' => 0.613,
	'14' => 0.851,
	'15' => 0.732,
	'16' => 0.828,
	'17' => 0.615,
	'18' => 0.804,
	'19' => 0.685,
	'20' => 0.584,
);	


if (defined $fasta_file)
{	print "$fasta_file \n";
	my $seq = readfasta($fasta_file);
	my @keys = keys %$seq;
	foreach my $k (@keys)
	{ 
		print "$k \n$seq->{$k} \n";
		
		my $crispr_seq = detectcrispr($seq->{$k});
		my $upper = scalar(@{$crispr_seq});
		for (my $i = 0; $i < $upper; $i++)
		{
			print $i;
#		print $i."\t".$crispr_seq->[$i]->{"sgRNAandPAM"}."\t".$crispr_seq->[$i]->{"sgRNA"}."\t".FirstScreen($crispr_seq->[$i]->{"sgRNA"})."\n";	
#		ScoreCalculate($crispr_seq->[$i]->{"sgRNA"});
#print ScoreCalculate($crispr_seq->[$i]->{"sgRNA"});
		}
	
	}

}			


####################################
########### Read Fasta  ############
sub readfasta 
{
   (my $file)=@_;
	my %sequence;
	my $header;
	my $temp_seq;
	
	open (IN, "<$file") or die "couldn't open the file $file $!";
	while (<IN>)
	{	
		chop;
		next if /^\s*$/;
		if ($_ =~ /^>/) # when see head line
		{
			if ($temp_seq)# if already has contents
			{
				$temp_seq="";
				$header= $_;
			}
			else # if meet for the first time
			{
				$header=$_;	
			}
		}
		else # when see the sequence line 
		{
		   s/\s+//g;
		   $temp_seq .= $_;
		   $sequence{$header}=$temp_seq; #update the contents
		}
	
	}
	
	return \%sequence
}
########### End Read Fasta  ########
####################################



####################################
########### Detect Crispr ##########
# Three Parameters: 
# 1. Sequence 
# 2. Start position (if from fasta file, probably position is not available)
# 3. sgRNA length (can be used from global variable $length) 
# Output: a reference to an array of hash

sub detectcrispr
{
	my ($seq,$chr,$chr_sta,$chr_end) = @_;
	my @targets;
	
	unless (defined $chr_sta)
	{
	$chr_sta = 0; # use 0 based system 
	$chr = 0;		
	}
	
	while ($seq =~ /(?=gg)/gi) #First process plus strand
	{  	
		my $gg_start= $chr_sta + pos($seq);
		
		my $crispr_sta = $chr_sta+pos($seq)-$length-1;
		my $crispr_end= $chr_sta +pos($seq)+1;
		
		my $fetch_start =  pos($seq)-$length-1;
		my $fetch_end   =  pos($seq)-1;
		
		if ($crispr_sta >= $chr_sta)
		{

			my $target=substr($seq,$fetch_start,$length);
			my $full_target =substr($seq,$fetch_start,$length+3);
			my $temp_target;
			$temp_target ->{"chr"} = $chr;
			$temp_target ->{"sta"} = $crispr_sta;
			$temp_target ->{"end"} = $crispr_end;
			$temp_target ->{"sgRNA"} = $target;
			$temp_target ->{"sgRNAandPAM"} = $full_target;
			$temp_target ->{"strand"} = "plus";
			
			next if GetCount($target,2) > 5;
			next if GetCount($target,0) > 1;
			print $target."\t";
			print GetCount($target,0)."\t";
			print GetCount($target,1)."\t";
			print GetCount($target,2)."\t";
			print GetCount($target,3)."\t";
			my($score,$flag)=ScoreCalculate($target);
			print "$score\t$flag\n";
			
			push @targets, $temp_target;	# Create Array of Hash
		} 		
	}	
	
	while ($seq =~ /(?=cc)/gi) #Next process the minus strand
	{  	
		my $gg_start= $chr_sta + pos($seq);
		
		my $crispr_sta = $chr_sta+pos($seq)+3;
		my $crispr_end= $chr_sta +pos($seq)+$length+2;
		
		my $fetch_start =  pos($seq)+3;
		my $fetch_end   =  pos($seq)+$length+2;
		
		if ($fetch_end <= length($seq))
		{

			my $target=complementary(substr($seq,$fetch_start,$length)); 			# Get the complementary sequence
			my $full_target =complementary(substr($seq,$fetch_start-3,$length+3));  # Get the complementary sequence
			my $temp_target;
			$temp_target ->{"chr"} = $chr;
			$temp_target ->{"sta"} = $crispr_sta;
			$temp_target ->{"end"} = $crispr_end;
			$temp_target ->{"sgRNA"} = $target;
			$temp_target ->{"sgRNAandPAM"} = $full_target;
			$temp_target ->{"strand"} = "negative";
			
			
			next if GetCount($target,2) > 5;
			next if GetCount($target,0) > 1;
			print $target."\t";
			print GetCount($target,0)."\t";
			print GetCount($target,1)."\t";
			print GetCount($target,2)."\t";
			print GetCount($target,3)."\t";
			my($score,$flag)=ScoreCalculate($target);
			print "$score\t$flag\n";
			push @targets, $temp_target;	
		} 		
	}	
	
	return \@targets;

}

######## End Detech Crispr #########
####################################



####################################
######  Start complementary ########

sub complementary {
  my ($sequence) = @_;
  my $reverse_sequence = reverse($sequence); 
  $reverse_sequence =~ tr/ACGTRYNX/TGCAYRNX/;
  return $reverse_sequence;
}

#######  End complementary #########
####################################



####################################
######## Start GetCount   ##########
sub GetCount
{
	my ($target_sequence,$mismatch) = @_;
	my ($target_sequence) = @_;
	my $url_base = 'http://gggenome.dbcls.jp/'.$species;
	my $perfect_url=$url_base.'/'.$mismatch.'/'.$target_sequence.'.json';
	my $perfect_ua = LWP::UserAgent->new();
	my $perfect_response = $perfect_ua->get($perfect_url);
	my $perfect_json = JSON->new;
	my $perfect_data = $perfect_json->decode($perfect_response->content());
	my $total_count=$perfect_data->{'summary'}[0]{'count'}+$perfect_data->{'summary'}[1]{'count'};
	return $total_count;
	
}
########  End GetCount #############
####################################





####################################
sub ScoreCalculate
{	
	my ($target_sequence) = @_;
	my $url_base = 'http://gggenome.dbcls.jp/'.$species;
	my $max_missmatch = "3"; # Max Mismatch is 3 - same as online tool
	my $url = $url_base.'/'.$max_missmatch.'/'.$target_sequence.'.json';
	my $ua = LWP::UserAgent->new();
	my $response = $ua->get($url);
	my $json = JSON->new;
	my $data = $json->decode($response->content());
	my @results = @{$data->{'results'}};
	my $totalscore =0;
	my $flag= 0;
	foreach my $result ( @results ) 
	{	
		my $snippet =$result->{"snippet"}; # Get the snippet 0 based? 1 based? 
		# First deal with positive strand
		if ($result->{"strand"} eq '+')
		{	
			my $start= $result->{"position"} - $result->{"snippet_pos"}; #the mismatch sequence
			my $length = $result->{"position_end"} - $result->{"position"} + 1;	
			
			my $pam_start=$result->{"position_end"}-$result->{"snippet_pos"}+1;
			
			my $return_sequence= substr($snippet,$start,$length);
			my $pam_sequence = substr($snippet,$pam_start,3);
			
#			print  $target_sequence. "\n";
#			print  $return_sequence. "\n";
#			print needlemanwunsch($target_sequence,$return_sequence)."\n";
			
			my @return=diffseq($target_sequence,$return_sequence);
#			print join " ", @return;
#			print "\n";
			next if scalar @return == 0;
			if (@return <3){
			if ($pam_sequence =~ /.GG/){$flag = 1};} 
			my $score = sgScore(@return);
			
#			print $score."\n";
			
			$totalscore += $score;
#			print $result->{"name"}.":".$result->{"position"}."-".$result->{"position_end"}."\n";
#			print $pam_sequence."\n";
#			print "\n";
		}elsif($result->{"strand"} eq '-')
		{			
			my $start= $result->{"position"} - $result->{"snippet_pos"}; #the mismatch sequence
			my $length = $result->{"position_end"} - $result->{"position"} + 1;	
			
			my $pam_start=$result->{"position"}-$result->{"snippet_pos"}-3;
			
			my $return_sequence= complementary(substr($snippet,$start,$length));
			my $pam_sequence = complementary(substr($snippet,$pam_start,3));
			
#			print  $target_sequence. "\n";
#			print  $return_sequence. "\n";
#			print needlemanwunsch($target_sequence,$return_sequence)."\n";
			my @return=diffseq($target_sequence,$return_sequence);
#			print join " ", @return;
#			print "\n";
			next if scalar @return == 0;
			if (@return <3){
			if ($pam_sequence =~ /.GG/){$flag = 1};} 	
			my $score = sgScore(@return);
			
#			print $score."\n";
			
			$totalscore += $score;
#			print $result->{"name"}.":".$result->{"position"}."-".$result->{"position_end"}."\n";
#			print $pam_sequence."\n";
#			print "\n";			
		}
		
				
	}
#	print "the final score is $totalscore \n";
	my $per_finalscore= (100/(100+$totalscore))*100;
#	print "the final score is $per_finalscore \n";	
#	print "the flag is $flag \n";
	return ($per_finalscore,$flag);
}


####################################
######## Start ParseUCSC   #########

sub ParseUCSC
{
	my ($chr,$sta,$end) = @_;
	my ($link) ="http://genome.cse.ucsc.edu:80/cgi-bin/das/mm10/dna?segment=$chr:$sta,$end";
	print "$link \n";
	my ($genefile) = get($link);
    my @DNA= grep {/^[acgt]*$/i;}split("\n",$genefile);
	my ($fasta_out)=$chr."_".$sta."_".$end.".fa";	
	return (@DNA);
}

########### End ParseUcsc  #########
####################################



####################################
########### Start diffseq ##########

sub diffseq 
{
	my ($align1, $align2) = needlemanwunsch(@_) ;
	my @mis;
	foreach (0..length($align2)-1)
	{
		if (substr($align1,$_,1) eq substr($align2,$_,1))
		{
			next;	
		}elsif(substr($align1,$_,1) eq '-') 
		{
			push (@mis,$_+1);
		}elsif(substr($align2,$_,1) eq '-')		
		{
			push (@mis,$_+1);
		}else
		{
			push (@mis,$_+1)
		}
}

return @mis;

}

########### End diffseq ############
####################################



####################################
########### Start sgScore ##########
# Receive an array as input 
sub sgScore
{	
	my $FirstTerm = 1;
	my $SecondTerm;
		
	my @array = @_;
	my $numberofarray = scalar @array;
	my $min = min @array;
	my $max = max @array;
	my $DTerm;
	if ($numberofarray ==1)
		{$DTerm=19;
		}
		else
		{$DTerm = ($max-$min)/($numberofarray-1)
		};
	my $ThirdTerm = 1/($numberofarray*$numberofarray);
	$SecondTerm = 1/((19-$DTerm)/19*4+1);

	foreach my $position (@array)
	{
		if (defined $table{$position})
		{
			$FirstTerm *= 1-$table{$position}
		}
		
	} 
	
	$FirstTerm = $FirstTerm* 100;
	my $Score=$FirstTerm*$SecondTerm*$ThirdTerm;
	
#	say "D: $DTerm";
#	say "1st $FirstTerm ";
#	say "2nd $SecondTerm";
#	say "3rd $ThirdTerm";
	
	$Score = sprintf("%.1f", $Score);
#	say "Score: $Score";
	return $Score;
}

####################################
###### NeedleManwunsch Algorithm ####

sub needlemanwunsch {

# Needleman-Wunsch Algorithm

# get sequences
my $seq1 = $_[0] // '' ;
my $seq2 = $_[1] // '' ;

# scoring scheme
my $MATCH    =  1; # +1 for letters that match
my $MISMATCH = -1; # -1 for letters that mismatch
my $GAP      = -1; # -1 for any gap

# initialization
my @matrix;
$matrix[0][0]{score}   = 0;
$matrix[0][0]{pointer} = "none";
for(my $j = 1; $j <= length($seq1); $j++) {
	$matrix[0][$j]{score}   = $GAP * $j;
	$matrix[0][$j]{pointer} = "left";
}
for (my $i = 1; $i <= length($seq2); $i++) {
	$matrix[$i][0]{score}   = $GAP * $i;
	$matrix[$i][0]{pointer} = "up";
}

# fill
for(my $i = 1; $i <= length($seq2); $i++) {
	for(my $j = 1; $j <= length($seq1); $j++) {
		my ($diagonal_score, $left_score, $up_score);
		
		# calculate match score
		my $letter1 = substr($seq1, $j-1, 1);
		my $letter2 = substr($seq2, $i-1, 1);		
		if ($letter1 eq $letter2) {
			$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
		}
		else {
			$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
		}
		
		# calculate gap scores
		$up_score   = $matrix[$i-1][$j]{score} + $GAP;
		$left_score = $matrix[$i][$j-1]{score} + $GAP;
		
		# choose best score
		if ($diagonal_score >= $up_score) {
			if ($diagonal_score >= $left_score) {
				$matrix[$i][$j]{score}   = $diagonal_score;
				$matrix[$i][$j]{pointer} = "diagonal";
			}
			else {
				$matrix[$i][$j]{score}   = $left_score;
				$matrix[$i][$j]{pointer} = "left";
			}
		} else {
			if ($up_score >= $left_score) {
				$matrix[$i][$j]{score}   = $up_score;
				$matrix[$i][$j]{pointer} = "up";
			}
			else {
				$matrix[$i][$j]{score}   = $left_score;
				$matrix[$i][$j]{pointer} = "left";
			}
		}
	}
}

my $align1 = "";
my $align2 = "";

my $j = length($seq1);
my $i = length($seq2);

while (1) {
	last if $matrix[$i][$j]{pointer} eq "none";
	
	if ($matrix[$i][$j]{pointer} eq "diagonal") {
		$align1 .= substr($seq1, $j-1, 1);
		$align2 .= substr($seq2, $i-1, 1);
		$i--; $j--;
	}
	elsif ($matrix[$i][$j]{pointer} eq "left") {
		$align1 .= substr($seq1, $j-1, 1);
		$align2 .= "-";
		$j--;
	}
	elsif ($matrix[$i][$j]{pointer} eq "up") {
		$align1 .= "-";
		$align2 .= substr($seq2, $i-1, 1);
		$i--;
	}	
}

$align1 = reverse $align1;
$align2 = reverse $align2;

return ($align1, $align2);
} ;


###### NeedleManwunsch Algorithm ####
####################################

