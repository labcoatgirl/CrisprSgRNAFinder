package MyModules::ReadBed;

use strict;
use lib '..';
use vars qw($VERSION);

use MyModules::FetchDNAfromUCSC qw/GetDNA_from_UCSC/;

$VERSION     = 1.00;

use base 'Exporter';
our @EXPORT_OK = qw/ReadBed/;


sub ReadBed
{
	my ($file,$genome) = @_;
	my @AoH;
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
		my $seq= GetDNA_from_UCSC($bedcolumns[0],$bedcolumns[1],$bedcolumns[2],$genome);
		push @AoH, { seq => $seq, chr => $bedcolumns[0], start => $bedcolumns[1], end => $bedcolumns[2]};
			
	}
	return \@AoH; # return array of hash 
}

1;