use strict;
use MyModules::FetchDNAfromUCSC qw/GetDNA_from_UCSC/;
use MyModules::ReadBed qw/ReadBed/;
use MyModules::ReadFasta qw/ReadFasta/;
use MyModules::DNAStuff qw/GetReverseComplementary  GetGCPercentage/;

use MyModules::ListAllsgRNA qw/ListAllsgRNA/;

use Data::Dumper;



use Term::ANSIColor;



# ----------------------------------------------------------------------
# Test ReadFasta Function
# ----------------------------------------------------------------------
my $fasta_file = "test_input_0.fasta";

my $seq_fasta = ReadFasta($fasta_file);
my $seq_num = scalar(@$seq_fasta);
print colored("Finish reading $fasta_file\n$seq_num sequence(s) have been read\n", "bold blue");

for my $item (@$seq_fasta) 
{
	print $$item{chr}."\t";
	print $$item{start}."\t";
	print $$item{end}."\t";
	print $$item{seq}."\n";
	
	print Dumper ListAllsgRNA($item);
}

	
# ----------------------------------------------------------------------
# Test ReadBed module and DNAStuff module
# ----------------------------------------------------------------------	
	
my $bed_file = "test_input_2.bed";
my $genome ="hg19";

my $seq_bed = ReadBed($bed_file,$genome);
my $seq_num = scalar(@$seq_bed);
print colored("Finish reading $bed_file\n$seq_num sequence(s) have been read\n", "bold blue");

for my $item (@$seq_bed) 
{
	print $$item{chr}."\t";
	print $$item{start}."\t";
	print $$item{end}."\t";
	print $$item{seq}."\n";
	my $sequence = $$item{seq};
	print GetReverseComplementary($sequence)."\t";
	print GetGCPercentage($sequence)."\n";
	
	
	
}

