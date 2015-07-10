package ReadBed;

use strict;
use Exporter;
use vars qw($VERSION);

$VERSION     = 1.00;


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

1: