use strict;
use MyModules::FetchDNAfromUCSC;
use MyModules::ReadBed;
use MyModules::ReadFasta qw/ReadFasta/;
use Term::ANSIColor;


my $fasta_file = "test_input_0.fasta";

my $seq = ReadFasta($fasta_file);
	my $seq_num = scalar(@$seq);
	print colored("Finish reading $fasta_file\n$seq_num sequence(s) have been read\n", "bold blue");
	
	for my $item (@$seq) 
	{
#		my $mytargets = detectsgRNA($item);
#		evaluatePAM($mytargets);
		print "=========================================\n$item\n";
	}

my $DNA = GetDNA_from_UCSC("chr3","8780000","8785000","mm10");
print $DNA;