package MyModules::QuerysgRNA;

use lib '..';
use strict;
use vars qw($VERSION);
use base 'Exporter';
use LWP::UserAgent;
use LWP::Simple;
use JSON;
use MyModules::DNAStuff qw/GetReverseComplementary  GetGCPercentage/;


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
return $total_count;

}