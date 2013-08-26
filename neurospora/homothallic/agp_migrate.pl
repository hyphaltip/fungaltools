#!/usr/bin/perl -w
use strict;
use File::Spec;
use Bio::DB::Fasta;

my $min_size = 100;
# currently we assume all velvet AGP files are all + stranded

my $topdir = '/shared/stajichlab/projects/neurospora_homothallic/assemblies';
my $velvet = File::Spec->catfile($topdir,'velvet/RunOn');
my $mercator = File::Spec->catfile($topdir,'mercator');


my %species;

opendir(V, $velvet) || die $!;
for my $file ( readdir(V) ) {
    next unless ($file =~ /(\S+)\.agp/);

    my $stem = $1;
    $stem =~ s/_woHC//;
    open(my $fh => File::Spec->catfile($velvet,$file)) || die $!;
    # these scf names refer to the NODE_XXXX names
    my %l;
    while(<$fh>) {
	next if /^#/;
	next if $l{$_}++;
	my ($scaf,@row) = split;
	my ($scaf_id);
	if( $scaf =~ /^scf_(\d+)/) {
	    $scaf_id = $1;
	} else {
	    die("cannot find parseable scaffold id in $scaf\n");
	}
	push @{$species{$stem}->{velvet}->{$scaf_id}}, [@row];
    }
}

opendir(M, $mercator) || die $!;
my %scaff_names;
for my $file ( readdir(M) ) {
    next unless ($file =~ /(\S+)\.agp/);
    my $stem = $1;
    open(my $ofh => ">$stem.EMBL.agp") || die $!;
    if( ! $species{$stem} ) { 
	warn("cannot find $stem in species set\n");
	next;
    }
    open(my $fh => File::Spec->catfile($mercator,$file)) || die "$file: $!";
    my %frag_order;	
    my @info;
    while(<$fh>) {
	next if /^#/;
	my ($scaffold,$scaf_start,$scaf_end, $order,
	    $type, 
	    @rest) = split;
	my $scaffold_name;
	my ($sn,$st);
	if($scaff_names{$scaffold} ) {
	    $scaffold_name = $scaff_names{$scaffold};
	    ($st,$sn) = split('_',$scaffold_name);	       
	} elsif( $scaffold =~ /^assembled(\d+)/) {
	    $scaff_names{$scaffold} = $scaffold_name = sprintf("SC_%d",$1);
	    $st = 'SC';
	    $sn = $1;
	} elsif( $scaffold =~ /^NODE_(\d+)/) {
	    $scaff_names{$scaffold} = $scaffold_name = sprintf("scf_%d",$1);
	    $st = 'scf';
	    $sn = $1;
	} else {
	    die("unmatched scaffold name $scaffold\n");
	}
	my $len = abs($scaf_end - $scaf_start);
	if( $len < $min_size ) {
	    push @info, [$st, $sn, 
			 $scaffold_name, $scaf_start,$scaf_end, $order, 
			 'N', $len, 'scaffold', 'yes'];
	    
	} else {
	    push @info, [$st, $sn, 
			 $scaffold_name, $scaf_start,$scaf_end, $order, $type,
			 @rest];
	}
    }
    open(my $outname => ">$stem.scaffold_name_lookup.dat") || die $!;
    for my $nm ( 
	# schwartzian transformation
	map { $_->[2] } 
	sort { $a->[0] cmp $b->[0] ||
			    $a->[1] cmp $b->[1] } 
		 map { [ split('_',$scaff_names{$_},2), $_] }
		 keys %scaff_names ) {
	print $outname join("\t", $scaff_names{$nm},$nm), "\n";
    }
    for my $d ( sort { $a->[0] cmp $b->[0] ||
			   $a->[1] <=> $b->[1] } @info ) {
	my ($st,$sn,
	    $scaffold_name, $scaf_start,$scaf_end, $order, $type,
	    @rest) = @$d;

	if( $type eq 'N' ) {
	    print $ofh join("\t", 
			    $scaffold_name,
			    $scaf_start, $scaf_end, 
			    ++$frag_order{$scaffold_name}, # insure '1' 1st time
			    'U',
			    100, # fixed size for 'U' type gaps
			    'scaffold',
			    'yes',
			    'align_xgenus'),"\n";
	} else {
	    my ($contig, $ctg_start, $ctg_end, $ctg_strand) = @rest;
	    my $ctg_id;
	    if( $contig =~ /NODE_(\d+)_/) {
		$ctg_id = $1;
	    } else {
		die("cannot parse ctg id from $contig\n");
	    }
	    if( ! exists $species{$stem}->{velvet}->{$ctg_id} ) {
		  die "no data for $ctg_id in velvet agp!\n";		  
	    }
	    my $rows = $species{$stem}->{velvet}->{$ctg_id};
	    my ($running_start, $running_end) = ($scaf_start, $scaf_end);
	    for my $row ( @$rows ) {
		if( $row->[3] eq 'N' ) {
		    print $ofh join("\t",
				    $scaffold_name,
				    $running_start, 
				    $running_start + $row->[4] -1,
				    ++$frag_order{$scaffold_name},
				    'N',
				    $row->[4],
				    'scaffold',
				    'yes',
				    'paired-ends'), "\n";		    
		    $running_start += $row->[4];
		} else {
		    print $ofh join("\t",
				    $scaffold_name,
				    $running_start,
				    $running_start + ($row->[6] -$row->[5]),
				    ++$frag_order{$scaffold_name},
				    'D',
				    $row->[4], $row->[5], $row->[6],
				    $row->[7]
			), "\n";
		    $running_start += ($row->[6] - $row->[5])+1;
		}
	    }
	}
    }
}



