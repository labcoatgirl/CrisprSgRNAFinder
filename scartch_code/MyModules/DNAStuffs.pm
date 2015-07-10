sub GetReverseComplementary {
  my ($sequence) = shift;
  my $reverse_sequence = reverse($sequence); 
  $reverse_sequence =~ tr/ACGTRYNX/TGCAYRNX/;
  return $reverse_sequence;
}


1;