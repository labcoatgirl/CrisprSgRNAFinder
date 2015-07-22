package MyModules::ReadFasta;

use strict;
use vars qw($VERSION);
use base 'Exporter';
use Term::ANSIColor;

$VERSION     = 1.00;

our @EXPORT_OK = qw/ReadFasta/;


sub ReadFasta 
{
    (my $file)= shift;
	my %sequence;
	my $header="";
	my $temp_seq;
	
	#suppose the fasta file contains multiple sequences;
	 
	open (IN, "<$file") or die "couldn't open the file $file $!";
	while (<IN>)
	{	
		$_ =~ s/\R//g;
		next if /^\s*$/; #skip empty line 
		if ($_ =~ /^>/)  #when see head line
		{	
			$header= $_;
			if ($sequence{$header}){print colored("#CAUTION: SAME FASTA HAS BEEN READ MULTIPLE TIMES.\n
			#CAUTION: PLEASE CHECK FASTA SEQUENCE:$header\n","red")};
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
	
	my @AoH;
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
		push @AoH,  { seq => $sequence{$k}, chr => $chr, start => $sta, end => $end};	
	}
	
	undef %sequence;
	
	return \@AoH;
}

1;