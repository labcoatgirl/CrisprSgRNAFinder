
use strict;
use lib '..';
use Data::Dumper;
use Term::ANSIColor;

use MyModules::FetchDNAfromUCSC qw/GetDNA_from_UCSC/;
use MyModules::ReadBed qw/ReadBed/;
use MyModules::ReadFasta qw/ReadFasta/;
use MyModules::DNAStuff qw/GetReverseComplementary  GetGCPercentage/;
use MyModules::ListAllsgRNA qw/ListAllsgRNA/;
use MyModules::QuerysgRNA qw/QuerysgRNA/;

# ----------------------------------------------------------------------
# Test ReadFasta Function
# ----------------------------------------------------------------------
my $genome ="hg19";

my $fasta_file = "test_input_1.fasta";

my $seq_fasta = ReadFasta($fasta_file);
my $seq_num = scalar(@$seq_fasta);
print colored("Finish reading $fasta_file\n$seq_num sequence(s) have been read\n", "bold blue");

# for my $item (@$seq_fasta)
# {
# 	print $item->{chr}."\t";
# 	print $item->{start}."\t";
# 	print $item->{end}."\t";
# 	print $item->{seq}."\n";
	
	
	# my @sgRNA_list = ListAllsgRNA($item,"minus");
# 	$item->{sglist}=  \@sgRNA_list;
#
# 	foreach (@sgRNA_list)
# 	{
# 		my $guide_seq = $_ ->{guide_seq};
# 		my $pam_seq   = $_ ->{PAMseq};
# 		my $query_out = QuerysgRNA($guide_seq,1,$genome);
# 		print "1 match number for $guide_seq $pam_seq is $query_out\n";
#
# 		my $query_out = QuerysgRNA($guide_seq,2,$genome);
# 		print "2 match number for $guide_seq $pam_seq is $query_out\n";
#
# 		my $query_out = QuerysgRNA($guide_seq,3,$genome);
# 		print "3 match number for $guide_seq $pam_seq is $query_out\n";
# 	}
#}

#foreach (@sgRNA_list) 
#{
#	foreach my $key (sort keys($_))
#	{
#	print $_ ->{$key}."\t";
#	}
#	print "\n";
#}


#print Dumper ($seq_fasta);
	
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
	print $item->{chr}."\t";
	print $item->{start}."\t";
	print $item->{end}."\t";
	print $item->{seq}."\n";
	my $sequence = $item->{seq};	
	
	my @sgRNA_list = ListAllsgRNA($item,"both");
	$item->{sglist}=  \@sgRNA_list;
	
	foreach (@sgRNA_list)
	{
		my $guide_seq = $_ ->{guide_seq};
		my $pam_seq   = $_ ->{PAMseq};
		my $query_out = QuerysgRNA($guide_seq,1,$genome);
		print "1 mismatch allowed: ";
		print colored($guide_seq,"bright_green");
		print colored($pam_seq,"red");
		print " $query_out\n";

		my $query_out = QuerysgRNA($guide_seq,2,$genome);
		print "2 mismatch allowed: ";
		print colored($guide_seq,"bright_green");
		print colored($pam_seq,"red");
		print " $query_out\n";
		
		my $query_out = QuerysgRNA($guide_seq,3,$genome);
		print "3 mismatch allowed: ";
		print colored($guide_seq,"bright_green");
		print colored($pam_seq,"red");
		print " $query_out\n";
		
		print "----------------------------------------------------------------------------------------------------------------------------------\n";
		
	}
	
}

