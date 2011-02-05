use warnings;
use strict;
my $srx2strain = 'SRX2strain.txt';
my $srx2srp = 'SRX2SRP.txt';

my %srx;
open(my $fh => $srx2strain) || die $!;
# put the SRX -> strain into the hash
while(<$fh>) {
    next if /^\s+$/;
    my (undef,$srx,@strain) = split;
    $srx{$srx}->{strain} = join(" ", @strain);
}
open($fh => $srx2srp) || die $!;
# now put in the SRP values using the same key (SRX) which was in 1st file
while(<$fh>) {
    my ($srx,$srp) = split;
    $srx{$srx}->{srp} = $srp;
}
# print out the joined group
for my $srx ( sort keys %srx ) {
    print join("\t", $srx{$srx}->{srp}, $srx{$srx}->{strain},$srx),"\n";
}
