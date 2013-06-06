use LWP::Simple;
use strict;
my $gi_list = '24475906,224465210,50978625,9507198';

#assemble the URL
my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $query='txid4751[Organism:exp]';
my $db = 'genome';
my $url = $base . "esearch.fcgi?db=$db&term=$query&rettype=acc&usehistory=y";

#post the URL
my $output = get($url);

#parse WebEnv, QueryKey and Count (# records retrieved)
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);
print "$output\n";
#retrieve data in batches of 500
open(OUT, ">fungalgenomes.dat");
my $retmax = 500;
for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
        my $efetch_url = $base ."efetch.fcgi?db=$db&WebEnv=$web";
        $efetch_url .= "&query_key=$key&retstart=$retstart";
        $efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
        my $efetch_out = get($efetch_url);
        print OUT "$efetch_out";
}
close OUT;

#my @ids;
#parse IDs retrieved
#while ($output =~ /<Id>(\d+?)<\/Id>/sg) {
#   push(@ids, $1);
#}

