#!/usr/bin/perl -w
use strict;
use DBI;
use Getopt::Long;
use Text::CSV;

my $dbfile = 'phenoImages.db';
my $create = ! ( -f $dbfile && ! -z $dbfile );
my $debug = 0;
my $infile;
my $interval = 25_000;
GetOptions(
	   'c|create!' => \$create,
	   'v|verbose!'=> \$debug,
	   'i|in|input:s'=> \$infile,
	   'ci|commitinterval:i' => \$interval,
	   );

if( ! $infile ) {
    $infile = shift @ARGV;
}

my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile",undef,undef, # no username/password
		       { AutoCommit => 0,
			 RaiseError => 1,
			 sqlite_see_if_its_a_number => 1} );
$dbh->do("PRAGMA foreign_keys = ON");

if( $create ) {
    warn("creating table\n") if ($debug);
    $dbh->do("DROP TABLE IF EXISTS phenotype_obs");
    $dbh->do("DROP TABLE IF EXISTS image");
    $dbh->do("DROP TABLE IF EXISTS locus");

    $dbh->do(
'CREATE TABLE IF NOT EXISTS locus (					      
	locus_id	INTEGER PRIMARY KEY AUTOINCREMENT,
        locus           varchar(24) NOT NULL UNIQUE,
        gene_name       varchar(24)
)');	

    $dbh->do('CREATE TABLE IF NOT EXISTS image (					      
	image_id	INTEGER PRIMARY KEY, -- previously defined filenaem we downloaded
	locus_id INTEGER NOT NULL,
	extension varchar(6) NOT NULL, -- jpg or png
	title		varchar(48) NOT NULL,
	media		varchar(24) NOT NULL,
	username  varchar(24),
	dateCreated	datetime,
	dateModified	datetime,
        FOREIGN KEY(locus_id) REFERENCES locus(locus_id)
				       				      
)');


    $dbh->do('CREATE TABLE IF NOT EXISTS phenotype_obs (					      
	po_id	INTEGER PRIMARY KEY AUTOINCREMENT,
        locus_id INTEGER NOT NULL,
        media           varchar(24) NOT NULL,
	ontologyTerm	varchar(48) NOT NULL,
	value		double(8) NOT NULL,
	dateCreated	datetime,
	dateModified	datetime,
        FOREIGN KEY(locus_id) REFERENCES locus(locus_id)
)');

    my $locus_ins_sth = $dbh->prepare('INSERT INTO locus (locus) VALUES (?)');
    
    my $locus_id_sth = $dbh->prepare('SELECT locus_id FROM locus where locus = ?');

    my $image_ins_sth = $dbh->prepare('INSERT INTO image (image_id, locus_id, extension, title, media,
username,dateCreated, dateModified) VALUES (?,?,?,?,?,?,?,?)');

    my $phenotype_ins_sth = $dbh->prepare('INSERT INTO phenotype_obs (locus_id, media, ontologyTerm, value, dateCreated, dateModified) VALUES (?,?,?,?,?,?)');

    my $csv = Text::CSV->new ()
	or die "Cannot use CSV: ".Text::CSV->error_diag ();
    open(my $fh => $infile) || die "cannot open $infile: $!\n";
    
    my $hdr = <$fh>;
    chomp($hdr);
    $csv->column_names (split(/,/,$hdr));
    my $i = 0;
    my %seen;


    while(my $row = $csv->getline_hr($fh)) {
	if( keys %$row != 17 ) {
	    warn("uneven row\n");
#	    warn(scalar @$row, " cols ", "and row is ",join(";",@$row),"\n");
	    last;
	}
	for my $d ( values %$row ) {
	    chomp($d);
	}
	if( $row->{'pi.title'} =~ s/\s+(plate|edge)(\s+\S+)?$/ $1/) {
	    # drop some extraneous info that prevents this from being a general assay field
	}
#	my $status = splice(@$row,6,1);
#	my $status2 = splice(@$row,10,1);
	my $locus = $row->{'al.locus'} || die "no locus found in the row\n";
	
	$locus_id_sth->execute($locus);
	my $r = $locus_id_sth->fetch;
	my $locus_id;
	if( $r && $r->[0] ) {
	    $locus_id = $r->[0];
	} else {
	    $locus_ins_sth->execute($locus);
	    $locus_id = $dbh->last_insert_id('','','locus','locus_id');
	    warn("$locus_id for $locus\n") if $debug;
	}
	my $fname = $row->{'pi.id'}. '.'.$row->{'pi.extension'};
	unless( $seen{$fname} ++ ) {
	    $image_ins_sth->execute($row->{'pi.id'},
				    $locus_id,
				    map { $row->{$_} }
				    qw(pi.extension pi.title p.name
				       pi.username pi.dateCreated pi.dateModified) );
	}
	unless( $seen{join(",",
			  $row->{'p.dateCreated'},
			   $row->{'po.ontologyTerm.text'},
			   $row->{'pm.value'},
			   $row->{'p.dateModified'} )}++ ) {
	    $phenotype_ins_sth->execute($locus_id,
					map { $row->{$_} } qw(p.name po.ontologyTerm.text pm.value p.dateCreated p.dateModified) );
	}
	$dbh->commit if( $i > 0 && $i % $interval == 0);
	last if $debug && $i > 1000;
	$i++;
    }
    $dbh->commit;
    my $counter_sth = $dbh->prepare("SELECT COUNT(*) FROM locus");
    $counter_sth->execute();
    my $counter_row = $counter_sth->fetch;
    warn( "INSERTED ", $counter_row->[0], " loci\n");
    $dbh->do("CREATE INDEX i_ontTerm ON phenotype_obs (ontologyTerm)");
}
 
my $query = $dbh->prepare("SELECT locus.locus, image.image_id, image.extension, image.title, image.media  FROM locus, image WHERE image.locus_id = locus.locus_id");
my ($locus,$image_id,$ext,$title,$media);
$query->execute;
$query->bind_columns(\$locus,\$image_id,\$ext,\$title,\$media);
while( $query->fetch ) {
    my $imgfile = sprintf("%s/%d.%s",$ext,$image_id,$ext);
    #print join("\t", $locus, $imgfile, $title, $media), "\n";
    if( ! -f $imgfile ) {
	warn("missing $imgfile\n");
    }
}

END {

    $dbh->disconnect;

}
