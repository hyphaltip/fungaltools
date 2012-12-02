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
    $dbh->do("DROP TABLE IF EXISTS phenotype");
    $dbh->do("DROP TABLE IF EXISTS ontology_term");
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
        fname           varchar(24) NOT NULL,
	title		varchar(48) NOT NULL,
	username  varchar(24),
	dateCreated	date,
        FOREIGN KEY(locus_id) REFERENCES locus(locus_id)
				       				      
)');
    $dbh->commit;
    $dbh->do('CREATE TABLE IF NOT EXISTS ontology_term (
         ontology_term_id   INTEGER PRIMARY KEY AUTOINCREMENT,
         classification varchar(24) NOT NULL,
         value varchar(48) NOT NULL
)');
    $dbh->do('CREATE INDEX i_onterm ON ontology_term (classification,value)');
    $dbh->commit;
    $dbh->do('CREATE TABLE IF NOT EXISTS phenotype (					      
	po_id	INTEGER PRIMARY KEY AUTOINCREMENT,
        locus_id INTEGER NOT NULL,
	ontology_term_id   INTEGER NOT NULL,
        obs_type varchar(24) NOT NULL,
        condition           varchar(24) NOT NULL,
        status          varchar(32) NOT NULL,	
        FOREIGN KEY(locus_id) REFERENCES locus(locus_id)
--        FOREIGN KEY(ontology_term_id) REFERENCES locus(ontology_term_id)
)');

    $dbh->commit;
    my $locus_ins_sth = $dbh->prepare('INSERT INTO locus (locus) VALUES (?)');
    
    my $locus_id_sth = $dbh->prepare('SELECT locus_id FROM locus where locus = ?');

    my $image_ins_sth = $dbh->prepare('INSERT INTO image (image_id, locus_id, extension, fname, title,username,dateCreated) VALUES (?,?,?,?,?,?,?)');

    my $phenotype_ins_sth = $dbh->prepare('INSERT INTO phenotype (locus_id, obs_type, condition, ontology_term_id, status) VALUES (?,?,?,?,?)');

    my $ont_ins_sth = $dbh->prepare('INSERT INTO ontology_term (classification,value) VALUES (?,?)');

    my $ont_id_sth = $dbh->prepare('SELECT ontology_term_id FROM ontology_term where classification = ? AND value = ?');

    my $csv = Text::CSV->new ({'sep_char' => "\t"})
	or die "Cannot use CSV: ".Text::CSV->error_diag ();
    open(my $fh => $infile) || die "cannot open $infile: $!\n";
    
    my $hdr = <$fh>;
    chomp($hdr);
    $csv->column_names (split(/\t/,$hdr));
    my $i = 0;
    my %seen;
    my %ont_terms;
    my %loci;
    while(my $row = $csv->getline_hr($fh)) {
	
	my $locus = $row->{'Locus'} || die "no locus found in the row\n";
	my $locus_id;
	unless( exists $loci{$locus} ) {
	    my $r = $locus_id_sth->fetch;	 
	    if( $r && $r->[0] ) {
		$locus_id = $r->[0];
	    } else {
		$locus_ins_sth->execute($locus);
		$locus_id = $dbh->last_insert_id('','','locus','locus_id');
		warn("$locus_id for $locus\n") if $debug;
	    }
	    $loci{$locus} = $locus_id;
	} else {
	    $locus_id = $loci{$locus};
	}
	my $fname = $row->{'Image File'};
	unless( $seen{$fname} ++ ) {
	    my ($image_id,$ext) = split(/\./,$fname);
	    $image_ins_sth->execute($image_id,
				    $locus_id,
				    $ext,
				    $fname,
				    $row->{'Image Title'},
				    $row->{'User'},
				    $row->{'Date'} );
	}
	
	my ($classif,@ont_values) =  split(':',$row->{'Ontology Term'});
	for my $val ( @ont_values ) {
	    my $ont_term_id;
	    if( ! exists $ont_terms{$classif}->{$val} ) {
		$ont_id_sth->execute($classif,$val);
		my $r = $ont_id_sth->fetch;
		if( $r && $r->[0] ) {
		    $ont_term_id = $ont_terms{$classif}->{$val} = $r->[0];
		} else {
		    $ont_ins_sth->execute($classif,$val);
		    $ont_term_id = $ont_terms{$classif}->{$val} = $dbh->last_insert_id('','','ontology_term','ontology_term_id');
		}
	    } else {
		$ont_term_id = $ont_terms{$classif}->{$val};
	    }
	    $phenotype_ins_sth->execute($locus_id,
					$row->{'Observation Type'},
					$row->{'Conditions'},
					$ont_term_id,
					$row->{'Status'});
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
}
 
my $query = $dbh->prepare("SELECT locus.locus, image.image_id, image.extension, image.title  FROM locus, image WHERE image.locus_id = locus.locus_id");
my ($locus,$image_id,$ext,$title);
$query->execute;
$query->bind_columns(\$locus,\$image_id,\$ext,\$title);
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
