package FetchDNAfromUCSC;

use strict;
use vars qw($VERSION);
#use LWP::UserAgent;
#use LWP::Simple;
use XML::Simple;

$VERSION     = 1.00;


sub Parse_UCSC
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