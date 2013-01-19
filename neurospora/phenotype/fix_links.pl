#!/usr/bin/perl -w
use strict;
use File::Spec;
my @dirs = qw(jpg png);

for my $dir ( @dirs ) {
    opendir(IN, $dir) || die $!;
    for my $file ( readdir(IN) ) {
	if( -l "$dir/$file" ) {
	    my $full_link = readlink("$dir/$file");
	    my (undef,undef,$lookup_name) = File::Spec->splitpath($full_link);
	    
	    my ($base,$ext) = split(/\./,$lookup_name);

	    my $target_file = $lookup_name;
	    if( ! -f $full_link || ! -f "$dir/$target_file" ) {
		print "$file $full_link $lookup_name\n";
		if( $ext ne $dir ) {
		    $target_file = "../$ext/$lookup_name";
		}
		unlink( "$dir/$file");
		my $success = symlink($target_file,"$dir/$file");
		if( ! $success ) {
		    warn "failed to symlink $target_file to $dir/$file\n";
		}
	    }
	}
    }
}
