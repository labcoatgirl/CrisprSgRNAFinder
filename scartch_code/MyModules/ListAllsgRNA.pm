package MyModules::ListAllsgRNA;

use lib '..';
use strict;
use vars qw($VERSION);
use base 'Exporter';
use MyModules::DNAStuff qw/GetReverseComplementary  GetGCPercentage/;


$VERSION     = 1.00;
our @EXPORT_OK = qw/ListAllsgRNA/;

sub ListAllsgRNA
{
	my $input = shift;
	my $seq  = $$input{seq};
	my $chr	 = $$input{chr};
	my $chr_sta = $$input{start};
	my $chr_end = $$input{end};
	
	my $lower_gc_threshold = 40.0;
	my $upper_gc_threshold = 60.0;
	
	my $guidesequence_length = 20; 
	
	my @sgRNAs_temp ; #initilize array of hash
	my @sgRNAs_final; 
	my %unique_check;
	
	
	while ($seq =~ /(?=GG)/gi) #First process + strand 
	{  	
		# An illustration for sgRNA position 
		# positions and sequences are as shown below:
		# if start position is not provided, use 0 based system
		# 0 - 1 - 2 - 3 - 4 - 5 - 6 - 7 
		# A - C - C - T - G - G - C - C
		# pos($seq) = 4 -> the G-G starts at position 4
		# suppose guide sequence lenght is 2 bp not including NGG (usually we have 20 nt for sequence before NGG) 
		# guideseq_sta should start at position 1
		# guideseq_sta = 4 - length(2) - 1   
		
		
		# In a different scenario where the start position is not 0
		# 1306 1307 1308 1309 1310 1311 1312 1313 1314 1315
		#  A    C     C    T   G     C    G    G    A    A
		#  0    1     2    3   4     5    6    7    8    9
		# pos($seq) = 6 -> the G-G starts at position 6
		# suppose guide sequence lenght is 2 bp not including NGG 
		# guideseq_sta should start at position 1309
		# guideseq_sta = 6 - length(2) - 1 + start = 6 - 2 - 1 + 1306 = 1309 
		# guideseq_sta = start + pos(seq) - length(sgRNA) - 1 
		# guideseq_end_with_PAM = start + pos(seq) + 1 
		
		
		
		# A general form is  
		# guideseq_sta = 0 + pos(seq) - length(sgRNA) - 1 
		# guideseq_end_with_PAM = 0 + pos(seq) + 1 
		
		my $guideseq_sta 			= $chr_sta + pos($seq) - $guidesequence_length - 1;
		my $guideseq_end   			= $chr_sta + pos($seq) - 2;
		my $guideseq_plus_PAM_end 	= $chr_sta + pos($seq) + 1;
		
		
       	# In evaluation of off-targets, we only use the guide sequence 
		# The whole sequnce (guide_sequence + PAM) is saved as guideseq_plus_PAM

		if ($guideseq_sta >= $chr_sta)
		{
			my $guideseq			=	substr($seq, $guideseq_sta-$chr_sta, $guidesequence_length);
			my $guideseq_plus_PAM 	=	substr($seq, $guideseq_sta-$chr_sta, $guidesequence_length+3);
			my $PAMseq 				=	substr($seq, $guideseq_end-$chr_sta+1, 3);
			
			my $temp_target;			
			$temp_target ->{"chr"} 					= $chr;
			$temp_target ->{"guide_seq_sta"} 		= $guideseq_sta;
			$temp_target ->{"guide_seq_end"} 		= $guideseq_end;
			$temp_target ->{"guide_seq"} 			= $guideseq;
			$temp_target ->{"guideseq_plus_PAM"} 	= $guideseq_plus_PAM;
			$temp_target ->{"PAMseq"}				= $PAMseq; 	
			$temp_target ->{"strand"} 				= "+";
			
			$unique_check{$guideseq} ++;
			
			push @sgRNAs_temp, $temp_target;	# Create a temporary Array of Hash
		} 		
	
	}	
	
	
	while ($seq =~ /(?=CC)/gi) #Next process the minus strand
	{  	
		# 12 - 13 - 14 -15-16- 17 -18 - 19
		# 0  - 1  - 2 - 3 - 4 - 5 - 6 - 7 
		# A  - C  - C - T - G - G - C - C
		# First CC start at position 1 - which is also the end position of GG
		
		my $guideseq_sta 			= $chr_sta + pos($seq) + $guidesequence_length + 2;
		my $guideseq_end   			= $chr_sta + pos($seq) + 3;
		my $guideseq_plus_PAM_end 	= $chr_sta + pos($seq) ;
		
		
		if ($guideseq_sta  <= $chr_sta + length($seq) -1 )
		{
			my $guideseq			=	substr($seq, $guideseq_end-$chr_sta , $guidesequence_length);
			my $guideseq_plus_PAM 	=	substr($seq, $guideseq_plus_PAM_end-$chr_sta, $guidesequence_length+3);
			my $PAMseq 				=	substr($seq, $guideseq_plus_PAM_end-$chr_sta, 3);
			
			$guideseq 				= 	GetReverseComplementary($guideseq);
			$guideseq_plus_PAM		=	GetReverseComplementary($guideseq_plus_PAM);
			$PAMseq					= 	GetReverseComplementary($PAMseq);
			
			my $temp_target;			
			$temp_target ->{"chr"} 					= $chr;
			$temp_target ->{"guide_seq_sta"} 		= $guideseq_sta;
			$temp_target ->{"guide_seq_end"} 		= $guideseq_end;
			$temp_target ->{"guide_seq"} 			= $guideseq;
			$temp_target ->{"guideseq_plus_PAM"} 	= $guideseq_plus_PAM;
			$temp_target ->{"PAMseq"}				= $PAMseq; 	
			$temp_target ->{"strand"} 				= "-";
			
			$unique_check{$guideseq} ++;
			
			push @sgRNAs_temp, $temp_target;	# Create a temporary Array of Hash
			
			# ----------------------------------------------------------------------
		} 	# End of if Loop
	}		# End of While Loop
	
	
	
	foreach (@sgRNAs_temp) 
	{
	my $current_seq = $_ ->{"guide_seq"}; 
	my $current_pam_seq = $_ ->{"PAMseq"};
	my $current_gc = GetGCPercentage($current_seq);
	print "current sequence is $current_seq $current_pam_seq $unique_check{$current_seq}  $current_gc\t";
	next if $unique_check{$current_seq} > 1;
	
	print "Bingo bigger \n" if $current_gc < $lower_gc_threshold;
	print "Bingo lesser \n" if $current_gc > $upper_gc_threshold;
	
	push @sgRNAs_final, $_;
	
	}
	
	undef(@sgRNAs_temp);
	
	return @sgRNAs_final;
}

