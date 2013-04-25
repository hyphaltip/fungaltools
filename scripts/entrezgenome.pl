#!/usr/bin/perl -w
use strict;
use Bio::DB::EUtilities;
use Bio::Tools::EUtilities::History;
use LWP::Simple;

my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $query='txid4751[Organism:exp]';
my $db = 'genome';
my $url = $base . "esearch.fcgi?db=$db&term=$query&rettype=acc&retmax=1000&usehistory=y";

my $output = get($url);

#parse WebEnv, QueryKey and Count (# records retrieved)
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);
my @ids;
while ($output =~ /<Id>(\d+?)<\/Id>/sg) {
   push(@ids, $1);
}
#my $term = 'txid4751[Organism:exp]';
my $verbose =1;
#my $eutil = Bio::DB::EUtilities->new(
#				     -verbose => $verbose,
#				     -eutil      => 'esearch',
#				     -term       => $term,
#				     -db         => 'genome',
#				     -retmax     => 10000,
#				     -email      => 'jason.stajich@gmail.com'); 

# please use your real email

  # eutil => any of esearch, esummary, elink
#my @ids = $eutil->get_ids(); # returns array or array ref of IDs
#print scalar @ids,"\n";
for my $id ( @ids ) {
    my $fetch = Bio::DB::EUtilities->new(
					 -verbose => $verbose,
					 -eutil   => 'efetch',
					 -term    => $id,
					 -db      => 'genome',
					 -retmode => 'xml',
					 -webenv  => $web,
					 -email   => 'jason.stajich@gmail.com');
    print $fetch->get_Response(),"\n";   
    last;
}
