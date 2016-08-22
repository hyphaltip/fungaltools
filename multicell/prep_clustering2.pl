#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(sum);

my $clusters = 10;
my $skip_shared = 1;
my %singlecell = map { $_ => 1 } qw(SCER SPOM SROS Bcer Bjap Csym Drad Dvul Ecol Mmaz Phor Reut Rsph Ssol);
my %multicell = map { $_ => 1 } qw(CCIN ANID NCRA HSAP ATHA);
my $multicellspct = scalar keys %multicell;

GetOptions(
    'skip!' => \$skip_shared,
    );

my $dir = shift || ".";
opendir(DIR,$dir) || die $!;
for my $file ( readdir(DIR) ) {
    my %dat;
    my %seen;    
    next unless ($file =~ /(\S+)\.gene_patterns\.tab$/);
    my ($sp) = ($1);
    open(my $in => "$dir/$file" ) || die $!;
    my $ofile = "$sp.gene_pattern.dat";
    my $ofile_notsingle = "$sp.gene_pattern_nosinglecell.dat";
    my $ofile_multicell = "$sp.gene_pattern_allmulticell.dat";
    open(my $out => ">$dir/$ofile") || die $!;
    open(my $out_s => ">$dir/$ofile_notsingle") || die $!;
    open(my $out_m => ">$dir/$ofile_multicell") || die $!;
    while(<$in>) {
	my ($gene, $pattern) = split;
	my @group = split(/,/,$pattern);
	for my $s ( @group ) {
	    $dat{$gene}->{$s}++;
	    $seen{$s}++;
	}
    }
    my @exp = sort keys %seen;
    print $out join(",",qw(GENE),@exp),"\n";
    print $out_s join(",",qw(GENE),@exp),"\n";
    print $out_m join(",",qw(GENE),@exp),"\n";
    for my $gene ( sort keys %dat ) {
	my $count = sum ( map { $dat{$gene}->{$_} || 0 } @exp);
	next if $skip_shared && $count == scalar @exp; # skip those found in everyone (to reduce num of rows)
	my $print_singlecell = 1;
	my $print_multicell  = $multicell{$sp} || 0;

	for my $s ( keys %{$dat{$gene}} ) {
	    if( $singlecell{$s} ) {
		$print_singlecell = 0;
		last;
	    }
	    if ( $multicell{$s} && $s ne $sp ) {
		$print_multicell++; # requesting s!=$sp because we don't want to double count for paralogs
	    }
	}
	my $row = join(",",$gene, map { $dat{$gene}->{$_} || 0 } @exp);
	if( $count ) {
	    print $out_s $row, "\n" if $print_singlecell;
	    print $out_m $row, "\n" if( $print_singlecell && 
					($print_multicell / $multicellspct) >= 0.70 );
	}
	print $out $row,"\n";
    }
    open(my $R => ">$dir/$sp.heatmap.R") || die $!;
    print $R <<EOF
library(gplots);
library(fastcluster);
EOF
    ;
    printf $R 'dat <-read.table("%s",header=T,sep=",",row.names=1)'."\n",$ofile;
    print $R <<EOF
gm <- data.matrix(dat);
# choose some nice blue colors
palette <- c('#ffffff','#0033BB')
h.rows <- hclust(dist(gm,method="manhattan"));
h.cols <- hclust(dist(t(gm),method="manhattan"));
EOF
    ;
    printf $R 'pdf("%s.matrix.pdf")'."\n", $sp;

if( 0 ) {
print $R <<EOF    
gene_heatmap <- heatmap(gm, Rowv=as.dendrogram(h.rows), Colv=as.dendrogram(h.cols),
col = palette, scale="none",margins=c(5,10));
EOF
    ;

    printf $R "mycl <- cutree(h.rows, k=%d);\n",$clusters;
print $R <<EOF
EOF
;


    for ( my $cluster_count = 1; $cluster_count < $clusters; $cluster_count++ ) {
	printf $R "clid <- c(%s);\n",$cluster_count;
print $R <<EOF
ysub <- gm[names(mycl[mycl%in%clid]),];
hrsub <- hclust(dist(ysub, method="manhattan"))                                        
# Select sub-cluster number (here: clid=c(1)) and generate corresponding dendrogram.
EOF
;

	printf $R "heatmap.2(ysub, main=\"cluster %d\",Rowv=as.dendrogram(hrsub),\n".
"Colv=as.dendrogram(h.cols), col=palette,\n".
'scale="none",trace="none",density.info="none");'."\n",
$cluster_count;	
    }
} # end comment out

    printf $R 'datsingl <-read.table("%s",header=T,sep=",",row.names=1)'."\n",$ofile_notsingle;
    print $R <<EOF
gm <- data.matrix(datsingl);
# choose some nice blue colors
h.rows <- hclust(dist(gm,method="manhattan"));
h.cols <- hclust(dist(t(gm),method="manhattan"));
EOF
    ;
    printf $R 'pdf("%s.matrix_notsinglecell.pdf")'."\n", $sp;

print $R <<EOF    
heatmap(gm, main="genes not in singlecell orgs",
	Rowv=as.dendrogram(h.rows), Colv=as.dendrogram(h.cols),
col = palette, scale="none",margins=c(5,10));
EOF
    ;


    printf $R 'datmult <-read.table("%s",header=T,sep=",",row.names=1)'."\n",$ofile_multicell;
    print $R <<EOF
gm <- data.matrix(datmult);
# choose some nice blue colors
h.rows <- hclust(dist(gm,method="manhattan"));
h.cols <- hclust(dist(t(gm),method="manhattan"));
EOF
    ;
    printf $R 'pdf("%s.matrix_inmulticell.pdf")'."\n", $sp;
    my $clusters_multi = 15;
    print $R "mycl <- cutree(h.rows, h=max(h.rows\$height)/2);\n";#,$clusters_multi;
    print $R "plot(h.rows); rect.hclust(h.rows, k=5, border=\"red\");\n";
    printf $R "mycolhc <- sample(rainbow(256));\n";
    print $R "mycolhc <- mycolhc[as.vector(mycl)];\n";
print $R <<EOF    
heatmap(gm, main="genes in multicellular orgs not in singlecell",
	Rowv=as.dendrogram(h.rows), Colv=as.dendrogram(h.cols),
	RowSideColors=mycolhc,
	col = palette, scale="none",margins=c(5,10));
EOF
    ;

if( 0 ) {
    for ( my $cluster_count = 1; $cluster_count <= $clusters_multi; $cluster_count++ ) {
	printf $R "clid <- c(%s);\n",$cluster_count;
print $R <<EOF
ysub <- gm[names(mycl[mycl%in%clid]),];
hrsub <- hclust(dist(ysub, method="manhattan"))                                        
# Select sub-cluster number (here: clid=c(1)) and generate corresponding dendrogram.
EOF
;
	printf $R "heatmap.2(ysub, main=\"cluster %d\",Rowv=as.dendrogram(hrsub),\n".
"Colv=as.dendrogram(h.cols), col=palette,\n".
'scale="none",trace="none",density.info="none",margins=c(5,10));'."\n",
$cluster_count;	
    }
}
    close($R);
}
