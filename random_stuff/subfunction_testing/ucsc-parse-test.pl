use LWP::Simple;


$genome="mm9";
sub Parse_UCSC
{
	my ($chr,$sta,$end) = @_;
	my ($link) ="http://genome.ucsc.edu/cgi-bin/das/"."$genome"."/dna?segment=$chr:$sta,$end";
	my $xml = XML::Simple->new();
	$data = $xml->XMLin(get($link));
	my $sequence = $data-> {SEQUENCE}->{DNA}->{content};
  	# chop $sequence;
	$sequence =~ s/\s//;
	return $sequence;
}

$a = Parse_UCSC(chr1,5000000,5000300);
 print $a;
 print $a;

