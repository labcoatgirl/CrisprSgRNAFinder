#!/usr/bin/perl

use strict;
use warnings;
use 5.006;

use Getopt::Long qw( GetOptionsFromString );



### --- constant values ---

my $L_PAM = 4; # length of PAM in search progress & output; 4 nt because of the exist of -NNGG
my $L = 12; # l_nonseed can be 7-9
my $L_PAM_INPUT = 3; # length of PAM in input target file; last 3 nt of each input sequence is PAM
my $NUM_notify = 200; # frequency to notify in scanning of the sequences in genome file (including scaffolds)

my @NT = qw( A C G T N ); # used to build the hash of PAM (2nd-4th nt)
my @nt = qw( a c g t n ); # used to build the hash of seed off-targets

# maxism mismatches (MM) allowed currently
my $MAX_opt_mm_seed = 6; # if > 4, it may be very slow or out of memory
my $DEFAULT_opt_mm_seed = 2;
my $MAX_opt_mm_nonseed = 255; # means no limit
my $MAX_opt_mm_all = 255; # means no limit
my @PAM_LEVELS = qw( A B C N );

# for paired mode
my $MAX_opt_distance = 1000;

# for target mode
my $MIN_LEN_SITE = 18; # minimum length of protospacer
my $MAX_LEN_SITE = 30; # maximum length of protospacer
my $MAX_LEN_TARGET = 1000;

my @LEN_SITE_PRIORITY; # which length of target site is prefered in target mode; generally, the smaller, the prior.
$LEN_SITE_PRIORITY[$_] = $_ foreach ( $MIN_LEN_SITE .. $MAX_LEN_SITE );
$LEN_SITE_PRIORITY[20] = -2; # 20-nt protospacer has 1st priority
$LEN_SITE_PRIORITY[19] = -1; # 19-nt protospacer has 2nt priority





### --- set default values, get and check options ---

my $opt_help;
my $opt_wait = $^O =~ /Win32/ ? 1 : 0; # wait before exit: -w=0, not wait; -w=1, wait.

my $opt_mode = 'single';
my $opt_output = 'csv';

my $opt_mm_seed;
my $opt_mm_nonseed;
my $opt_mm_all;
my $opt_pam_level = 'A';

# for paired mode
my $opt_distance = 100; # paired mode: maximum distance allowed

# for target mode
my $opt_require5g = 'yes'; # target mode: whether the 5'-G (for T7 RNA polymerase) is required
my $opt_length = '19-20'; # target mode: length of target to search

# options of file names
my ( $opt_target, @opt_genome, $opt_exon );


my $errmsg = "Use '$0 -h' for quick help; for more information, please see README or visit http://eendb.zfgenetics.org/casot/";

my $helpmsg = qq{
===== CasOT: CRISPR/Cas9 system (Cas9/gRNA) Off-Targeter =====

--- Usage ---
	$0 -t=<target_site_file> -g=<genome_sequence_file> [<options>]
e.g.,
	$0 -t=example -g=GRCh37 -e=r73h.gtf -o=tab -s=1 -n=4 -p=B

For more information, please see README or visit http://eendb.zfgenetics.org/casot/

--- Options ---

	-h: show this help

	-m=<searching_mode>: single, paired, or target (default: $opt_mode)
	-t=<target_site_file>
	-g=<genome_sequence_file>
	-e=<exon_annotation_file>
	-o=<output_format>: csv or tab (defualt: $opt_output)

	-s=<maximum_number_of_mismatches_allowed_in_the_seed_region>: 0-$MAX_opt_mm_seed (nt) (default: $DEFAULT_opt_mm_seed; >4 is not suggested)
	-n=<maximum_number_of_mismatches_allowed_in_the_non-seed_region>: 0-$MAX_opt_mm_nonseed (nt) (default: $MAX_opt_mm_nonseed)
	-p=<level_of_PAM_type_allowed>: A, B, C or N (default: $opt_pam_level)
		Level A: NGG
		Level B: NGG + NAG
		Level C: NGG + NAG + NNGG
		Level N: (no limit)

	-d=<maximum_distance_of_paired_off-target>: 0-$MAX_opt_distance (default: $opt_distance)

	-r=<require_5'-G_or_not_in_candidate_target>: yes or no (default: $opt_require5g)
	-l=<length_range_allowed_for_candidate_target>: two number between 18 and 30 (default: $opt_length)

};



my %options= (
	'help' => \$opt_help,
	'wait=i' => \$opt_wait,

	'mode=s' => \$opt_mode,
	'output=s' => \$opt_output,

	'target=s' => \$opt_target,
	'genome=s' => \@opt_genome,
	'exon=s' => \$opt_exon,

	'seed=i' => \$opt_mm_seed,
	'nonseed=i' => \$opt_mm_nonseed,
	'all=i' => \$opt_mm_all,
	'pam=s' => \$opt_pam_level,

	'distance=i' => \$opt_distance,

	'require5g=s' => \$opt_require5g,
	'length=s' => \$opt_length,
);

if ( @ARGV ) {
	GetOptions( %options );
}
else {
	print
		"\n" .
		"Please input the option string (in one line), and press <Enter>.\n" .
		"For quick help, input '-h' and press <Enter>.\n" .
		"\n" .
		"The web-based opt-generator can be used to generate option string; if you are using Microsoft Windows, copy the string from generator, right-click the title (NOT the body) of the command prompt windows, then choose `Edit' - `Paste'.\n" .
		"\n";
	my $optstring;
	while ( !$optstring ) {
		print "# Option string: ";
		$optstring = <STDIN>;
		$optstring =~ s{^\s+|\s+$}{}g;
	}
	(undef, undef) = GetOptionsFromString( $optstring, %options );
}


print "\n";

err( $helpmsg ) if $opt_help;

$opt_mode = lc $opt_mode;
err( "ERROR - Unknown mode: -m $opt_mode (supported mode: single, paired, target).\n$errmsg" )
	unless grep /^$opt_mode$/, qw( single paired target );

$opt_output = lc $opt_output;
err( "ERROR - Unknown output format: -o $opt_output (supported mode: csv, tab ).\n$errmsg" )
	unless grep /^$opt_output$/, qw( csv tab );

if ( defined $opt_mm_all ) {
	err( "ERROR - Cannot use -a & -s/-n at the same time.\n$errmsg" )
		if defined($opt_mm_seed) || defined($opt_mm_nonseed); # cannot use -a & -s/-n at the same time
	err( "ERROR - Wrong option of number of mismatches allowed in the whole protospacer: -a $opt_mm_all (supported: 0-$MAX_opt_mm_all).\n$errmsg" )
		if $opt_mm_all < 0 || $opt_mm_all > $MAX_opt_mm_all;
	$opt_mm_seed = $opt_mm_nonseed = $opt_mm_all;
}
else { # if -a is not specified
	$opt_mm_seed = $DEFAULT_opt_mm_seed if !defined $opt_mm_seed;
	$opt_mm_nonseed = $MAX_opt_mm_nonseed if !defined $opt_mm_nonseed;
	$opt_mm_all = $MAX_opt_mm_all;
	err( "ERROR - Wrong option of number of mismatches allowed in the seed region: -s $opt_mm_seed (supported: 0-$MAX_opt_mm_seed).\n$errmsg" )
		if $opt_mm_seed < 0 || $opt_mm_seed > $MAX_opt_mm_seed;
	print "WARNING - Option '-s 5' or '-s 6' need very large memeory and long running time.\n\n"
		if $opt_mm_seed > 4;
	err( "ERROR - Wrong option of number of mismatches allowed in the non-seed region: -n $opt_mm_nonseed (supported: 0-$MAX_opt_mm_nonseed).\n$errmsg" )
		if $opt_mm_nonseed < 0 || $opt_mm_nonseed > $MAX_opt_mm_nonseed;
}

$opt_pam_level = uc $opt_pam_level;
err( "ERROR - Wrong option of PAM type allowed: -p $opt_pam_level (suported: " . ( join ", ", @PAM_LEVELS ) . ").\n$errmsg" )
	unless grep /^$opt_pam_level$/, @PAM_LEVELS;

err( "ERROR - Wrong option of distance allowed of two potential off-targets in the paired mode: -d $opt_distance (supported: 0-$MAX_opt_distance).\n$errmsg" )
	if $opt_distance < 0 || $opt_distance > $MAX_opt_distance;

$opt_require5g =
	$opt_require5g =~ /^y/i	?	1	:
	$opt_require5g =~ /^n/i	?	0	: err( "ERROR - Wrong option of require 5'-G: -r $opt_require5g (supported: yes or no).\n$errmsg" )
;
$opt_length =~ /^\s*(\d+)-(\d+)\s*$/;
err( "ERROR - Wrong option of target length region allowed: -l $opt_length (usage: -l <min>-<max>; example: -l 19-20; <min> and <max> should between $MIN_LEN_SITE & $MAX_LEN_SITE).\n$errmsg" )
	if !$1 || !$2 || $1<$MIN_LEN_SITE || $1>$MAX_LEN_SITE || $2<$MIN_LEN_SITE || $2>$MAX_LEN_SITE;
my ( $min_len_t, $max_len_t ) = $1<$2 ? ($1, $2) : ($2, $1);


err( "ERROR - No target file provided.\n$errmsg" ) if !$opt_target;
err( "ERROR - No genome sequence file provided.\n$errmsg" ) if !@opt_genome;





### --- open input & output files ---

open my $fseq, '<', $opt_target or err( "ERROR - Cannot open target file: $opt_target" );
open my $fexon, '<', $opt_exon or err( "ERROR - Cannot open annotation file: $opt_exon" )
	if $opt_exon;
foreach ( @opt_genome ) {
	open my $fgenome, '<', $_ or err( "ERROR - Cannot open genome file: $_" );
	close $fgenome;
}

( my $tfilename = $opt_target ) =~ s{\.fa$}{}; # remove the path and the suffix .fa for output directory name
my @gfilename = @opt_genome;
$_ =~ s{^.*/|^.*\\|\.fa$}{}g foreach @gfilename; # remove the path and the suffix .fa for output directory name
my $outpath = "$tfilename-" . (join '+', @gfilename) . "-s$opt_mm_seed";
$outpath .= "n$opt_mm_nonseed" if $opt_mm_nonseed < $MAX_opt_mm_nonseed;
$outpath .= "a$opt_mm_all" if $opt_mm_all < $MAX_opt_mm_all;
$outpath .= "p$opt_pam_level" if $opt_pam_level gt 'A';
mkdir $outpath, 0755 unless -e $outpath;
# err( "ERROR - Cannot make output directory: $outpath\n$errmsg" ) if $!; # always true in Windows?

open my $fstat, '>', "$outpath/_stat.txt" or err( "ERROR - Cannot open output statistic file: $outpath/_stat.txt" );
open my $fsite, '>', "$outpath/_sites.txt" or err( "ERROR - Cannot open output target file of the target mode: $outpath/_sites.txt" )
	if $opt_mode eq 'target';

# output files (%fout) will be opened in an individual section





### --- define global variables used in multiple subs ---

# created in section of reading target file & used in all subs
my ( @target_names, @target_pair_names ); # for sorting
my ( %targets, %l_nonseed, %l_all );

# created in section of building the seed hash & non-seed arrays & used in scan()
my %seed_ot; # candidate off-targets of seed regions of all targets
my %nts_of_nonseed; # nucleotide arrays of non-seed regions of each targets
# created in section of building the seed hash & non-seed arrays & used in scan() & outputpair()
my $longest_nonseed = 0;
my $mm_nonseed_pre_zero = ''; # whether a zero should be added before count of mismatches in non-seed region; A23 vs. A203

# created in scan() & used in findpair()
# should be cleared before scanning new chr
my ( %half_loc, %half_info );

# created in section of reading annotation file & used in getexon()
my %exon;
# used in getexon() and should be cleared before calling of scan() or outputpair()
my ( $prev_chr, $exon_cnt_of_chr, $idx_exon );

# created in findpair() & used in outputpair()
my @otpair;

# created in scan() & used in the section of outputting statistic information
my ( %stat, %stat_pair ) ;

# seperator in output
my $S = $opt_output eq 'tab' ? "\t" : ",";





### --- section of reading target file ---

# mode: target
if ( $opt_mode eq 'target' ) {
	no warnings 'numeric';
	my ( $seq );
	while ( <$fseq> ) {
		chomp;
		next if /^$/ || /^>/;
		$_ = uc;
		s/[^ACGTN]//g;
		$seq .= $_;
	}
	my $len_seq = length $seq;
	err( "ERROR - the sequence is too long (> $MAX_LEN_TARGET nt).\n$errmsg" ) if $len_seq > $MAX_LEN_TARGET;
	my $regx = $opt_require5g
		?	"G.{$min_len_t,$max_len_t}GG"
		:	".{" . ($min_len_t+1) . "," . ($max_len_t+1) . "}GG";
	while ( $seq =~ /$regx/g ) {
		my $len_site = length $&;
		my $name = pos($seq) . 'f';
		$targets{$name} = $&
			if  !$targets{$name}  ||  $LEN_SITE_PRIORITY[ $len_site - $L_PAM_INPUT ] < $LEN_SITE_PRIORITY[ length($targets{$name}) - $L_PAM_INPUT ];
		pos($seq) -= $len_site-1;
	}
	( $seq = reverse $seq ) =~ tr/ACGT/TGCA/;
	while ( $seq =~ /$regx/g ) {
		my $len_site = length $&;
		my $name = ( $len_seq - pos($seq) +1 ) . 'r';
		$targets{$name} = $&
			if  !$targets{$name}  ||  $LEN_SITE_PRIORITY[ $len_site - $L_PAM_INPUT ] < $LEN_SITE_PRIORITY[ length($targets{$name}) - $L_PAM_INPUT ];
		pos($seq) -= $len_site-1;
	}
	@target_names = sort { $a <=> $b } ( keys %targets );
	print $fsite ">$_\n$targets{$_}\n" foreach ( @target_names );
	close $fsite;
	use warnings 'numeric';
}

# mode: single 
elsif ( $opt_mode eq 'single' ) {
	my $unnamed_cnt = 0;
	my ( $tname, $tseq ) = ( '', '' );
	while ( <$fseq> ) {
		chomp;
		next if /^$/;
		if ( /^>/ ) {
			( $tname = substr $_, 1 ) =~ s/\s.*//;
			err( "ERROR - More than one sequences have the same name: $tname" ) if $targets{$tname};
			push @target_names, $tname;
		}
		else {
			err( "ERROR - Sequence without name: $_\n$errmsg" ) if !$tname;
			$tseq = checkseq($_);
			err( "ERROR - Sequence with wrong format: >$tname - $_" ) if !$tseq;
			$targets{$tname} = $tseq;
			$tname = '';
		}
	}
}

# mode: paired
elsif ( $opt_mode eq 'paired' ) {
	my ( %target_names_F, %target_names_R );
	my ( $tname, $tseq ) = ( '', '' );
	my $idx = 1; # for sorting the pairs in input order
	while ( <$fseq> ) {
		chomp;
		next if /^$/;
		if ( /^>/ ) {
			( $tname = substr $_, 1 ) =~ s/\s.*//;
			err( "ERROR - Sequence names in paired mode should end with '_#F' or '_#R': $tname\n$errmsg" ) if $tname !~ /_#(F|R)$/;
			err( "ERROR - More than one sequences have the same name: $tname" ) if $targets{$tname};
			$target_names_F{ substr $tname, 0, -3 } = $idx++ if $tname =~ /_#F$/;
			$target_names_R{ substr $tname, 0, -3 } = 1 if $tname =~ /_#R$/;
		}
		else {
			err( "ERROR - Sequence without name: $_\n$errmsg" ) if !$tname;
			$tseq = checkseq($_);
			err( "ERROR - Wrong target sequence: >$tname - $_\n$errmsg" ) if !$tseq;
			$targets{$tname} = $tseq;
			$tname = '';
		}
	}
	foreach ( keys %target_names_F ) {
		err( "ERROR - Unpaired sequence: >${_}_#F - " . $targets{"${_}_#F"} . "\n$errmsg" )
			if !$target_names_R{$_};
	}
	foreach ( keys %target_names_R ) {
		err( "ERROR - Unpaired sequence: >${_}_#R - " . $targets{"${_}_#R"} . "\n$errmsg" )
			if !$target_names_F{$_};
	}
	@target_pair_names = sort { $target_names_F{$a} <=> $target_names_F{$b} } keys %target_names_F;
	@target_names = map( ( "${_}_#F", "${_}_#R" ), @target_pair_names );
}

close $fseq;



sub checkseq {
	my( $s ) = @_;
	$s = uc $s;
	$s =~ s/\s+//g;
	return $s if
		$s =~ /^[ACGTN]*$/
		&& length $s >= $MIN_LEN_SITE + $L_PAM_INPUT
		&& length $s <= $MAX_LEN_SITE + $L_PAM_INPUT
		&& $s =~ /GG$/
	;
	return 0;
}





### --- opening output files ---

my %fout;

my $ext = $opt_output eq 'tab' ? 'txt' : 'csv';

foreach ( @target_names ) {
	open $fout{$_}, '>', "$outpath/$_.$ext" or err( "ERROR - Cannot open output file: $outpath/$_.$ext" );
	my $fout = $fout{$_};
	print $fout join $S, ( '# Location', qw( Site Target Mm.Type PAM Mm.All ) );
	print $fout "${S}Exon_info" if $opt_exon;
	print $fout "\n";
}

foreach ( @target_pair_names ) {
	open $fout{$_}, '>', "$outpath/$_.$ext" or err( "ERROR - Cannot open pairing output file: $outpath/$_.$ext" );
	my $fout = $fout{$_};
	print $fout join $S, ( '# Location', qw( Whole_site Target_pair Distance Location_+ Site_+ Target_+ Mm.Type_+ Location_- Site_- Target_- Mm.Type_- ) );
	print $fout "${S}Exon_info" if $opt_exon;
	print $fout "\n";
}
 




### --- outputting ready message ---

my $readymsg =
	"Summary of options\n" .
	"------------------\n" .
	"Input:\n" .
	"\tSearching mode: $opt_mode\n" .
	"\tTarget " . ($opt_mode eq 'target' ? "sequence" : "site") . " file: $opt_target\n" .
	"\tGenome file(s): " . (join ", ", @opt_genome) . "\n" .
	($opt_exon ? "\tExon annotation file: $opt_exon\n" : '' )
;
$readymsg .=	
	"Mismatches:\n" .
	"\tMaximum number of mismatches allowed in the seed region: $opt_mm_seed nt\n" .
	"\tMaximum number of mismatches allowed in the non-seed region: " . ($opt_mm_nonseed==$MAX_opt_mm_nonseed ? "no limit" : "$opt_mm_nonseed nt") . "\n" .
#	($opt_mm_all==$MAX_opt_mm_all ? "" : "\tMaximum number of mismatches allowed in the whole protospacer: $opt_mm_all\n") .
	"\tAllowed level of PAM type: " . ($opt_pam_level eq 'N' ? "no limit" : "Level $opt_pam_level") . "\n"
;
$readymsg .=
	"For `paired' mode:\n" .
	"\tMaximum distance allowed between the two potential off-target sequences: $opt_distance nt\n"
		if $opt_mode eq 'paired'
;
$readymsg .=
	"For `target' mode:\n" .
	"\tGuanine in the first position is required? $opt_require5g\n" .
	"\tAllowed range of protospacer length of candidate sites: $opt_length nt\n"
		if $opt_mode eq 'target'
;
$readymsg .=
	"Output:\n" .
	"\tOutput format: $opt_output\n" .
	"\tOutput path: $outpath\n" .
	"\tThe statistic file: $outpath/_stat\n"
;
$readymsg .= "\tThe file of candidate target sites: $outpath/_sites\n"
	if $opt_mode eq 'target';
$readymsg .=
	"\n\n" .
	"Progress\n" .
	"--------\n" .
	"\n"
;
print $readymsg;

my $stime = time();
print "# Start time: " . localtime() . "\n\n";

sub outtime {
	my $t = shift;
	my $s = '';
	$s .= $t>=3600	?	sprintf "%2dh ", $t/3600	:	'    ';
	$t %= 3600;
	$s .= $t>=60	?	sprintf "%2dm ", $t/60		:	'    ';
	$t %= 60;
	$s .= sprintf "%2ds", $t;
	return $s;
}





### --- section of reading annotation file and building a hash of exon information ---

if ( $opt_exon ) {
	print outtime(time()-$stime) . "\t# Read annotation file ...\n";
	my ( $gid, $gname );
	while ( <$fexon> ) {
		next unless /^[^\t]+\t[^\t]+\texon\t/;
		my ( $chr, undef, $feature, $start, $end, undef, undef, undef, $info ) = split "\t";
		if ( $info =~ /gene_id "([^"]+)".*gene_name "([^"]+)"/ ) {
			$gid = $1;
			$gname = $2;
		}
		elsif ( $info =~ /gene_id "([^"]+)"/ ) {
			$gid = $1;
			$gname = '';
		}
		else {
			$gid = $gname = '';
		}
		push @{ $exon{$chr} }, { start => $start, end => $end, gid => $gid, gname => $gname };
	}
	$exon{$_} = [   sort   { $a->{'start'} <=> $b->{'start'} }   @{ $exon{$_} }   ]
		foreach ( keys %exon );
	close $fexon;
}





### --- section of building the hash of seed & arrays of non-seed

# build the hash of potential off-targets of the seed regions of all target sequences,
# build the arrays of nucleotides of non-seed region,
# and record the lengthes of non-seed regions and the whole sequences of all targets.

print outtime(time()-$stime) . "\t# Building off-target sequences pool ...\n";

my %seed_ot_f; # temporary hash, the key will be mixed with nucleotides in upper- or lowercase

foreach ( @target_names ) {
	$l_all{$_} = length $targets{$_};
	$l_nonseed{$_} = $l_all{$_} - $L_PAM_INPUT - $L;
	my $seed = substr $targets{$_}, $l_nonseed{$_}, $L;
	change_nt( $_, $seed, $opt_mm_seed, 1 );
	my $nonseed = substr $targets{$_}, 0, $l_nonseed{$_};
	$longest_nonseed = length($nonseed) if length($nonseed) > $longest_nonseed;
	$nts_of_nonseed{$_} = [ split '', $nonseed ];
}
$mm_nonseed_pre_zero = '0' if $longest_nonseed>=10 && $opt_mm_nonseed>=10 && $opt_mm_all-$opt_mm_seed>=10;

sub change_nt {
	my( $k, $s, $num_to_change, $start_pos ) = @_;
	$seed_ot_f{$s} .= "###$k\t+\t" . ($opt_mm_seed - $num_to_change) . "\t$s";
	if ( $num_to_change ) {
		for my $loc ( $start_pos .. $L ) {
			my $left = substr $s, 0, $loc-1;
			my $mid = lc substr $s, $loc-1, 1;
			my $right = substr $s, $loc;
			for my $nt ( @nt ) {
				next if $nt eq $mid;
				my $new_s = "$left$nt$right";
				change_nt( $k, $new_s, $num_to_change-1, $loc+1 );
			}
		}
	}
	return;
}

foreach ( keys %seed_ot_f ) {
	my $k = uc;
	$seed_ot{$k} .= $seed_ot_f{$_};

	# reverse complementary sequence of the OT will also be searched
	( $k = reverse $k ) =~ tr/ACGT/TGCA/;
	( my $rev_seq = reverse $_ ) =~ tr/ACGTacgt/TGCAtgca/;
	( my $v = $seed_ot_f{$_} ) =~ s{\t\+\t}{\t-\t}g;
	$v =~ s{$_}{$rev_seq}g;
	$seed_ot{$k} .= $v;
}

undef %seed_ot_f; # for saving RAM



# --- build hash of all types of PAM (the 2nd - 4th nt) ---

my ( %pam_desc, %pam_level );
foreach my $nt1 ( @NT ) {
	foreach my $nt2 ( @NT ) {
		foreach my $nt3 ( @NT ) {
			$pam_level{"${nt1}${nt2}${nt3}"} = 'N';
			$pam_desc{"${nt1}${nt2}${nt3}"} = '-';
		}
	}
}
foreach my $nt ( @NT ) {
	$pam_level{"${nt}GG"} = 'C';
	$pam_desc{"${nt}GG"} = 'NNGG';
}
foreach my $nt ( @NT ) {
	$pam_level{"AG$nt"} = 'B';
	$pam_desc{"AG$nt"} = 'NAG';
	$pam_level{"GG$nt"} = 'A';
	$pam_desc{"GG$nt"} = 'NGG';
}





### --- main(): read and treat genome files ---

print outtime(time()-$stime) . "\t# Searching off-targets in genome ...\n";

$prev_chr = '';

my $seqcount;

# NOTE: it is supposed that there is no number or spacer in sequence lines of the genome file, and all lines in this file ended with \n but not \r\n
foreach my $current_genome ( @opt_genome ) {
	open my $fgenome, '<', $current_genome or err( "ERROR - Cannot open genome file: $current_genome" );
	print "\n" . outtime(time()-$stime) . "\t# Begin to search $current_genome\n";
	my( $seqname, $seq );
	$seqcount = 0;
	OUTER:
	while ( <$fgenome> ) {
		chomp;
		next OUTER if ( /^$/ );
		if ( /^>/ ) {
			$seqname = seqname( $_ );
			$seq = '';
			INNER:
			while (1) {
				if ( eof $fgenome ) { # the last sequence should also be scanned
					scan( $seqname, \$seq );
					last OUTER;
				}
				$_ = <$fgenome>;
				if ( /^>/ ) {
					scan( $seqname, \$seq );
					redo OUTER;
				}
				else {
					chomp;
					$seq .= uc $_;
				}
			}
			# end of INNER
		}
	}
	print outtime(time()-$stime) . "\t# All $seqcount sequences in $current_genome searched.\n";
	close $fgenome;
}


print "\n" . outtime(time()-$stime) . "\t# Outputing statistic information ...\n\n";
outputstat();

print "# End time: " . localtime() . " (" . outtime(time()-$stime) . ")\n\n";

close $fout{$_} foreach ( @target_names );
close $fout{$_} foreach ( @target_pair_names );
close $fstat;


print "\tOutput path: $outpath\n\n";

if ( $opt_wait ) {
	print "Press <Enter> to continue ...\n\n";
	<>;
}

1;
##### END of the program ######

sub seqname {
	my ( $name ) = @_;
	$name =~ s/^>([^ ]*) ?.*/$1/;
	print outtime(time()-$stime) . "\t  Searching chromosome: $name ...\n" if $name =~ /^[1-9]|^chr|^X|^Y/i; # chromosome, not scaffold
	print outtime(time()-$stime) . "\t  Searched $seqcount sequences.\n" if $seqcount && ( $seqcount % $NUM_notify == 0 );
	$seqcount ++;
	return $name;
}



sub err {
	my ( $msg ) = @_;
	print "\n$msg\n\n";
	if ( $opt_wait ) {
		print "Press <Enter> to continue ...\n";
		<>;
	}
	die "\n";
}





##### subs #####

### --- scanning genome ---

sub scan {
	my( $seqname, $r_seq ) = @_;
	%half_loc = %half_info = (); # the hash of location should be empty for each new chromosome

	my $l_seq = length $$r_seq;
	our ( $site, $pos ); # using 'our' is quicker than using 'my' here
	open my $fh, '<', $r_seq; # seek & read is quicker than substr
	for my $pos ( 0 .. $l_seq-$L-1 ) {
		seek $fh, $pos, 0;
		read $fh, $site, $L;

		if ( $seed_ot{$site} ) {
			my @tags = split '###', $seed_ot{$site};
			foreach ( @tags ) {
				next if !$_; # $_ may be ''

				my ( $tname, $strand, $mm_seed_cnt, $site ) = split "\t", $_;
				my $site_pam;

				# check PAM level first
				if ( $strand eq '+' ) {
					$site_pam = substr $$r_seq, $pos+$L+1, $L_PAM-1;
				}
				else {
					($site_pam = reverse substr $$r_seq, $pos-$L_PAM, $L_PAM-1) =~ tr/ACGT/TGCA/;
				}
				my $pam_level = $pam_level{$site_pam} or next; # no enough nts in start or end of the sequence for PAM
				next if $pam_level gt $opt_pam_level;
				my $pam_desc = $pam_desc{$site_pam};

				# if PAM allowed, get non-seed region and the whole PAM, and calculate the location for output
				my ( $site_nonseed, $tstart, $tend );
				if ( $strand eq '+' && $pos>=$l_nonseed{$tname} ) {
					$site_nonseed = substr $$r_seq, $pos-$l_nonseed{$tname}, $l_nonseed{$tname};
					$site_pam = substr $$r_seq, $pos+$L, $L_PAM;
					$tstart = $pos - $l_nonseed{$tname} + 1;
					$tend = $pos + $L + $L_PAM;
				}
				elsif ( $strand eq '-' && $pos<=$l_seq-$L-$l_nonseed{$tname} ) {
					($site_nonseed = reverse substr $$r_seq, $pos+$L, $l_nonseed{$tname}) =~ tr/ACGT/TGCA/;
					($site_pam = reverse substr $$r_seq, $pos-$L_PAM, $L_PAM) =~ tr/ACGT/TGCA/;
					$tstart = $pos - $L_PAM + 1;
					$tend = $pos + $L + $l_nonseed{$tname};
					($site = reverse $site) =~ tr/ACGTacgt/TGCAtgca/;
				}
				else {
					next; # no enough nts in start or end of the sequence for non-seed region
				}

				# calculate the mismatches in non-seed region
				my @nts_of_site_nonseed = split '', $site_nonseed;
				my $mm_nonseed_cnt =  0;
				$site_nonseed = '';
				for ( 0 .. $l_nonseed{$tname}-1 ) {
					if ( $nts_of_nonseed{$tname}[$_] ne $nts_of_site_nonseed[$_] ) {
						$mm_nonseed_cnt ++;
						$site_nonseed .= lc $nts_of_site_nonseed[$_];
					}
					else {
						$site_nonseed .= $nts_of_site_nonseed[$_];
					}
				}

				if ( $mm_nonseed_cnt <= $opt_mm_nonseed && $mm_seed_cnt + $mm_nonseed_cnt <= $opt_mm_all ) {
					my $whole_site = "${site_nonseed}_$site-$site_pam";
					my $mm_type = $mm_nonseed_cnt<10
						?	"$pam_level$mm_seed_cnt$mm_nonseed_pre_zero$mm_nonseed_cnt"
						:	"$pam_level$mm_seed_cnt$mm_nonseed_cnt";
					my $fout = $fout{$tname};
					print $fout join $S, (
						"$seqname:$tstart-$tend:$strand",
						"$whole_site",
						"$tname",
						"$mm_type",
						"$pam_level:$pam_desc",
						($mm_seed_cnt + $mm_nonseed_cnt)
					);
					print $fout getexon( $seqname, $tstart, $tend )
						if $opt_exon;
					print $fout "\n";
					$stat{$tname}{$mm_type} ++;

					if ( $opt_mode eq 'paired' ) {
						my ( $pairname, $halfname ) = split '_#', $tname;
						if ( $strand eq '+' ) {
							$pos += $L + $L_PAM -1;
						}
						else {
							$pos -= $L_PAM;
						}
						push @{ $half_loc{$pairname}{"$halfname$strand"} }, $pos;
						push @{ $half_info{$pairname}{"$halfname$strand"} }, "$tstart\t$tend\t$whole_site\t$mm_type";
					}
				} # end of output the off-target line

			}
		} # end of the existence test of $seed_ot{$site}

	}

	if ( $opt_mode eq 'paired' ) {
		@otpair = (); # clear the paired results for each chr.
		findpair( $seqname, $fh, 'F+', 'R-', 'Heter.' );
		findpair( $seqname, $fh, 'R+', 'F-', 'Heter.' );
		findpair( $seqname, $fh, 'F+', 'F-', 'Homo.' );
		findpair( $seqname, $fh, 'R+', 'R-', 'Homo.' );
		$prev_chr = ''; # getexon() has been used for single off-targets; it should be cleared before re-used for paired off-targets
		outputpair( $seqname );
	}

	close $fh;

} # end of the scan()





sub findpair {
	my ( $seqname, $fh, $locgroup_fwd, $locgroup_rev, $dimer_type ) = @_; # $dimer_type is not used currently; it may be used for TalenOT

	## NOTE: F or R means the forward and reverse single target site in a input pair
	## _fwd or _rev means the single off-target site located in the forward or reverse strand in a found potential off-target group

	my( $half_fwd, $strand_fwd ) = split '', $locgroup_fwd;
	my( $half_rev, $strand_rev ) = split '', $locgroup_rev;

	foreach my $t ( @target_pair_names ) {

		my $loc_fwd = $half_loc{$t}{$locgroup_fwd};
		my $loc_rev = $half_loc{$t}{$locgroup_rev};
		my $MIN = -( $opt_distance + $l_all{"${t}_#$half_fwd"} + $l_all{"${t}_#$half_rev"} );
		my $MAX = $opt_distance;

		next if ! defined $loc_fwd->[0] || ! defined $loc_rev->[0];
		my $idx_fwd = 0;
		my $idx_rev = 0;
		OUTER:
		for my $idx_fwd ( 0 .. $#{$loc_fwd} ) {
			my $pos_fwd = $loc_fwd->[$idx_fwd];
			while ( $loc_rev->[$idx_rev] < $pos_fwd+$MIN+1 ) {
				$idx_rev ++;
				last OUTER if !$loc_rev->[$idx_rev];
			}
			my $tmp_idx_rev = $idx_rev;
			while ( $loc_rev->[$tmp_idx_rev]  && $loc_rev->[$tmp_idx_rev] <= $pos_fwd+$MAX+1 ) {

				# output this pair

				my ( $start_fwd, $end_fwd, $whole_site_fwd, $mm_type_fwd ) = split "\t", $half_info{$t}{$locgroup_fwd}[$idx_fwd];
				my ( $start_rev, $end_rev, $whole_site_rev, $mm_type_rev ) = split "\t", $half_info{$t}{$locgroup_rev}[$tmp_idx_rev];

				( my $whole_site_rev_rev = reverse $whole_site_rev ) =~ tr/ACGTacgt/TGCAtgca/;

				my ( $start, $end, $whole_site, $distance );
				$start = $start_fwd < $start_rev ? $start_fwd : $start_rev;
				$end = $end_fwd > $end_rev ? $end_fwd : $end_rev;

				if ( $end_fwd < $start_rev ) {
					$distance = $start_rev - $end_fwd - 1;
					$whole_site = "(${t}_#$half_fwd $strand_fwd) ${whole_site_fwd}..N(" . $distance . ")..$whole_site_rev_rev (${t}_#$half_rev $strand_rev)";
				}
				elsif ( $start_fwd > $end_rev ) {
					$distance = $start_fwd - $end_rev - 1;
					$whole_site = "(${t}_#$half_rev $strand_rev) ${whole_site_rev_rev}..N(" . $distance . ")..$whole_site_fwd (${t}_#$half_fwd $strand_fwd)";
				}
				else { # overlapped
					seek $fh, $start-1, 0;
					read $fh, $whole_site, $end-$start+1;
					$distance = length($whole_site_fwd) + length($whole_site_rev) - ($end-$start+1) - 4; # each of the two whole sites has a '_' and a '-'
					$whole_site =
						$start_fwd < $start_rev
						?	"(${t}_#$half_fwd $strand_fwd) $whole_site (${t}_#$half_rev $strand_rev) (overlapped $distance)"
						:	"(${t}_#$half_rev $strand_rev) $whole_site (${t}_#$half_fwd $strand_fwd) (overlapped $distance)"
					;
					$distance = -$distance;
				}

				push @otpair, {
					t => $t,
					start => $start,
					end => $end,
					info => join $S, (
						"$seqname:$start-$end",
						"$whole_site",
						"$t",
						"$distance",
						"$seqname:$start_fwd-$end_fwd:$strand_fwd",
						"$whole_site_fwd",
						"${t}_#$half_fwd",
						"$mm_type_fwd",
						"$seqname:$start_rev-$end_rev:$strand_rev",
						"$whole_site_rev",
						"${t}_#$half_rev",
						"$mm_type_rev"
					)
				};

				$tmp_idx_rev ++;
			} # end of outputting the off-target line

		}
	} # end of each input target pair in this chromosome

} # end of findpair()





### --- output pairing off-targets ---

sub outputpair {
	my ( $seqname ) = @_;
	@otpair = sort { $a->{'start'} <=> $b->{'start'} || $a->{'end'} <=> $b->{'end'} } @otpair;
	foreach ( @otpair ) {
		my $fout = $fout{ $_->{'t'} };
		print $fout $_->{'info'};
		print $fout getexon( $seqname, $_->{'start'}, $_->{'end'} )
			if $opt_exon;
		print $fout "\n";
	}
}



sub getexon {
	my ( $chr, $tstart, $tend ) = @_;
	if ( !$prev_chr || $prev_chr ne $chr ) { # $prev_chr is a global variable, as well as $idx_exon & $exon_cnt_of_chr
		$idx_exon = 0;
		$prev_chr = $chr;
		$exon_cnt_of_chr = defined $exon{$chr} ? $#{ $exon{$chr} } : -1;
	}
	while ( $idx_exon <= $exon_cnt_of_chr ) {
		my $exon = $exon{$chr}[$idx_exon];
		# target at the left of current exon
		if ( $tend < $exon->{'start'} ) {
			last;
		}
		# target overlapped with current exon
		elsif ( $tstart <= $exon->{'end'} ) {
			return "$S$exon->{gname} ($exon->{gid})";
		}
		# target at the right of current exon
		else {
			$idx_exon ++;
		}
	}
	return "\t";
}







### --- output statistic data of single target ---
# there is no statistic data for paired target

sub outputstat {

	no warnings qw(uninitialized);

	foreach my $tname ( @target_names ) {
		print $fstat "\n--- $tname ---\n";
		if ( ! $stat{$tname} ) {
			print $fstat "no potential off-targets\n";
			next;
		}
		foreach ( sort keys %{ $stat{$tname} } ) {
			print $fstat "$_\t$stat{$tname}{$_}\n";
		}
	}

	print $fstat "\n\n\n------ Summary ------\n\n";
	print $fstat "Type\t", ( join "\t", @target_names ), "\n";
	foreach my $pam_level ( @PAM_LEVELS ) {
		last if $pam_level gt $opt_pam_level;
		print $fstat "\n";
		foreach my $mm_seed ( 0 .. $opt_mm_seed ) {
			print $fstat "\n";
			my $tmp_opt_mm_nonseed = $opt_mm_all - $mm_seed;
			$tmp_opt_mm_nonseed = $opt_mm_nonseed if $opt_mm_nonseed < $tmp_opt_mm_nonseed;
			$tmp_opt_mm_nonseed = $longest_nonseed if $longest_nonseed < $tmp_opt_mm_nonseed;
			foreach my $mm_nonseed ( 0 .. $tmp_opt_mm_nonseed ) {
				my $mm_type = $mm_nonseed<10
					?	"$pam_level$mm_seed$mm_nonseed_pre_zero$mm_nonseed"
					:	"$pam_level$mm_seed$mm_nonseed";
				print $fstat "$mm_type\t", ( join "\t", map { $stat{$_}{$mm_type} || '' } @target_names ), "\n";
			}
		}
	}

}
