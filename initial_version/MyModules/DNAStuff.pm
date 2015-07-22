package MyModules::DNAStuff;

use strict;
use vars qw($VERSION);
use base 'Exporter';

$VERSION     = 1.00;

our @EXPORT_OK = qw/GetReverseComplementary GetGCPercentage/;



sub GetReverseComplementary {
  my ($sequence) = shift;
  my $reverse_sequence = reverse($sequence); 
  $reverse_sequence =~ tr/ACGTRYNX/TGCAYRNX/;
  return $reverse_sequence;
}

sub GetGCPercentage {
	my ($sequence) = shift;
    my $sequencelength = length $sequence;
    my $GCcount = $sequence =~ tr/GC//;
    my $GCpercentage = ($GCcount+0.0)/($sequencelength+0.0);
	return $GCpercentage;
}


1;