package MyModules::QuerysgRNA;

use lib '..';
use strict;
use vars qw($VERSION);
use base 'Exporter';
use LWP::UserAgent;
use LWP::Simple;
use JSON;
use 5.010;
use MyModules::DNAStuff qw/GetReverseComplementary  GetGCPercentage/;
use MyModules::DiffCompare qw/DiffCompare/;

$VERSION     = 1.00;
our @EXPORT_OK = qw/QuerysgRNA/;

sub QuerysgRNA{

my ($target_sequence,$mismatch,$genome) = @_;
my $url_base = 'http://gggenome.dbcls.jp/'.$genome;
my $query_url;

if ($mismatch eq 0){
	$query_url=$url_base.'/'.$target_sequence.'.json'
}else{
	$query_url=$url_base.'/'.$mismatch.'/'.$target_sequence.'.json';
}

my $query_useragent = LWP::UserAgent->new();
my $query_response = $query_useragent->get($query_url);
my $query_json = JSON->new;
my $query_data = $query_json->decode( $query_response->content() );

my $total_count=$query_data->{'summary'}[0]{'count'}+$query_data->{'summary'}[1]{'count'};

my $results = $query_data->{'results'};


my $true_count;
my $high_risk_count;

foreach my $item( @$results) { 
	
    if ($item->{strand} eq "+"){
		# say "I am plus";
		my $guide_seq_sta =  $item->{position} - $item->{snippet_pos}; 
		my $guide_seq_end =  $item->{position_end} - $item->{snippet_pos};
		my $guide_seq_len =  $item->{position_end} - $item->{position} + 1;
		
		my $guideseq			=	substr($item->{snippet}, $guide_seq_sta, $guide_seq_len);
		my $pam					=   substr($item->{snippet}, $guide_seq_end+1, 3);
		
		# say "$guideseq $pam";
		
		if ($pam =~/.GG/){
			$true_count ++;
			my @array = DiffCompare($target_sequence,$guideseq);
			my $riskflag = 1; 
			foreach (@array){
				if ($_ > 8){$riskflag= undef}
			}
			if($riskflag){$high_risk_count ++}
			
		}
		
    }elsif ($item->{strand} eq "-"){ 
		# say "I am minus";
		my $guide_seq_sta =  $item->{position} - $item->{snippet_pos}; 
		my $guide_seq_end =  $item->{position_end} - $item->{snippet_pos};
		my $guide_seq_len =  $item->{position_end} - $item->{position} + 1;
		
		my $guideseq			=	GetReverseComplementary (substr($item->{snippet}, $guide_seq_sta, $guide_seq_len ));
		my $pam					=   GetReverseComplementary (substr($item->{snippet}, $guide_seq_sta-3, 3));
		# say "$guideseq $pam";
		if ($pam =~/.GG/){
			$true_count ++;
			my @array = DiffCompare($target_sequence,$guideseq);
			my $riskflag = 1; 
			foreach (@array){
				if ($_ > 8){$riskflag= undef}
			}
			if($riskflag){$high_risk_count ++}	
		}
    }
	
	
	
    # fields are in $item->{Year}, $item->{Quarter}, etc.
}

$total_count = $total_count-1;
$true_count  = $true_count-1;
if ($high_risk_count >=1){
	$high_risk_count = $high_risk_count-1;
}else{
	$high_risk_count =0
	
}

# say "true count is $true_count --over";
return "total mismatch count is $total_count | followed by NGG is $true_count| highrisk $high_risk_count";

#Try to add a more precise function to test whether the sequence follow a NGG


}