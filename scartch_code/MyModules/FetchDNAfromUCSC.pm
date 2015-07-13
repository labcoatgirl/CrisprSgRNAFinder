package MyModules::FetchDNAfromUCSC;

use strict;
use vars qw($VERSION);
use base 'Exporter';
use LWP::Simple qw/get/;
use XML::Simple;

$VERSION = 1.00;

our @EXPORT_OK = qw/GetDNA_from_UCSC/;

sub GetDNA_from_UCSC
{
	my ($chr,$sta,$end,$genome) = @_;
	my ($link) ="http://genome.ucsc.edu/cgi-bin/das/"."$genome"."/dna?segment=$chr:$sta,$end";
	my $xml = XML::Simple->new();
	my $data = $xml->XMLin(get($link));
	my $sequence = $data-> {SEQUENCE}->{DNA}->{content};
	$sequence =~ s/(\s+)//gi;
	$sequence = uc($sequence);
	return $sequence;
}

1;