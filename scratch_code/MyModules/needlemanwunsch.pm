
my ($seq1, $seq2) = @ARGV;
use Data::Dumper;

sub diff_compare{
my $seq1 = $_[0]; # The first sequence is expected to be the sgRNA 
my $seq2 = $_[1]; # The second sequence is expected to be the one with missmatch to designed sgRNA


# ----------------------------------------------------------------------
# Needleman-Wunsch Algorithm
#
# Code from: Korf I, Yandell M, Bedell J. BLAST. O'Reilly (2003)
# p.44 Example 3-1. Trace-back with Needleman-Wunsch algorithm
# Download Example Code: http://examples.oreilly.com/9780596002992/
# ----------------------------------------------------------------------


# scoring scheme
my $MATCH    =  1; # +1 for letters that match
my $MISMATCH = -1; # -1 for letters that mismatch
my $GAP      = -1; # -1 for any gap

# initialization
my @matrix;
$matrix[0][0]{score}   = 0;
$matrix[0][0]{pointer} = "none";
for(my $j = 1; $j <= length($seq1); $j++) {
	$matrix[0][$j]{score}   = $GAP * $j;
	$matrix[0][$j]{pointer} = "left";
}
for (my $i = 1; $i <= length($seq2); $i++) {
	$matrix[$i][0]{score}   = $GAP * $i;
	$matrix[$i][0]{pointer} = "up";
}

# fill
for(my $i = 1; $i <= length($seq2); $i++) {
	for(my $j = 1; $j <= length($seq1); $j++) {
		my ($diagonal_score, $left_score, $up_score);
		
		# calculate match score
		my $letter1 = substr($seq1, $j-1, 1);
		my $letter2 = substr($seq2, $i-1, 1);		
		if ($letter1 eq $letter2) {
			$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
		}
		else {
			$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
		}
		
		# calculate gap scores
		$up_score   = $matrix[$i-1][$j]{score} + $GAP;
		$left_score = $matrix[$i][$j-1]{score} + $GAP;
		
		# choose best score
		if ($diagonal_score >= $up_score) {
			if ($diagonal_score >= $left_score) {
				$matrix[$i][$j]{score}   = $diagonal_score;
				$matrix[$i][$j]{pointer} = "diagonal";
			}
			else {
				$matrix[$i][$j]{score}   = $left_score;
				$matrix[$i][$j]{pointer} = "left";
			}
		} else {
			if ($up_score >= $left_score) {
				$matrix[$i][$j]{score}   = $up_score;
				$matrix[$i][$j]{pointer} = "up";
			}
			else {
				$matrix[$i][$j]{score}   = $left_score;
				$matrix[$i][$j]{pointer} = "left";
			}
		}
	}
}

my $align1 = "";
my $align2 = "";

my $j = length($seq1);
my $i = length($seq2);

while (1) {
	last if $matrix[$i][$j]{pointer} eq "none";
	
	if ($matrix[$i][$j]{pointer} eq "diagonal") {
		$align1 .= substr($seq1, $j-1, 1);
		$align2 .= substr($seq2, $i-1, 1);
		$i--; $j--;
	}
	elsif ($matrix[$i][$j]{pointer} eq "left") {
		$align1 .= substr($seq1, $j-1, 1);
		$align2 .= "-";
		$j--;
	}
	elsif ($matrix[$i][$j]{pointer} eq "up") {
		$align1 .= "-";
		$align2 .= substr($seq2, $i-1, 1);
		$i--;
	}	
}

$align1 = reverse $align1;
$align2 = reverse $align2;
print "$align1\n";
print "$align2\n";

my @missmatch_positions ;
my @missmatch_type;
my $base_position =1;

foreach (0..length($align2)-1)
{
	if (substr($align1,$_,1) eq substr($align2,$_,1) )
	{	$base_position++;
		next; 
	}elsif(substr($align1,$_,1) eq '-')
	{	push @missmatch_positions,$base_position;
		push @missmatch_type,"insertion";
	}elsif(substr($align2,$_,1) eq '-')
	{	push @missmatch_positions,$base_position;
		push @missmatch_type,"deletion";
		$base_position++;
	}else{
		push @missmatch_positions,$base_position;
		push @missmatch_type,"mismatch";
		$base_position++;
	}

}

print Dumper (@missmatch_positions)."\n";
print Dumper (@missmatch_type)."\n";

}

diff_compare($seq1, $seq2);

