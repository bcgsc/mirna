#!/usr/bin/env perl
use strict;
use Getopt::Std;
use Pod::Usage;
use File::Basename;
use DBI;

use vars qw($opt_m $opt_o $opt_g $opt_v);
getopts("m:o:g:v:");

my $usage = "$0 -m mirbase_db -o species_code -g genome_version -v mirbase version name\n";
die "$usage" unless $opt_m && $opt_o && $opt_g && $opt_v;

my $db = $opt_m;
my $species = $opt_o;
my ($dbname, $dbhost, $dbuser, $dbpass) = get_db($db);
my $dbh_mirbase = DBI->connect("DBI:mysql:database=$dbname;host=$dbhost", $dbuser, $dbpass, {AutoCommit => 0, PrintError => 1}) || die "Could not connect to database: $DBI::errstr";

#get mirbase species code from organism code
my $species_code = $dbh_mirbase->selectrow_array("SELECT auto_id FROM mirna_species WHERE organism = '$species'");
die "$species_code organism code not found in database" unless defined $species_code;

#print header
print "miRNA_ID\tmiRBase_Ver\tAccession\tGenomic_Coords\tPrecursor_Seq\tMature_Coords\tMature_Accession\tAlt_Mature_Coords\tAlt_Mature_Accession\tStar_Coords\tStar_Accession\n";

#array of arrays, each subarray's indices correspond to:
#mirna id, mirna accession, coordinates, sequence, mature coordinates, mature MIMAT, alternate mature coordinates, alt mature MIMAT, star coordinates, star MIMAT
my @mirnas;
my ($mirna, $accession, $chr, $start, $end, $strand, $seq, $mature, $mature_start, $mature_end, $mature_acc);
my $sth = $dbh_mirbase->prepare("
	SELECT m.mirna_id, m.mirna_acc, co.xsome, co.contig_start, co.contig_end, co.strand, m.sequence, ma.mature_name, mature_from, mature_to, ma.mature_acc
FROM
mirna m
JOIN mirna_chromosome_build co ON m.auto_mirna = co.auto_mirna
JOIN mirna_pre_mature p ON m.auto_mirna = p.auto_mirna
JOIN mirna_mature ma ON p.auto_mature = ma.auto_mature
WHERE m.auto_species = '$species_code'
ORDER BY m.mirna_id, co.xsome, co.contig_start, co.contig_end, ma.mature_name
");
$sth->execute();
$sth->bind_columns(\$mirna, \$accession, \$chr, \$start, \$end, \$strand, \$seq, \$mature, \$mature_start, \$mature_end, \$mature_acc);
my $prev_mirna = "";
my $prev_coord = "";
while ($sth->fetchrow_arrayref()) {
	my $coord = "$opt_g:$chr:$start-$end:$strand";
	my $mature_coord = "$mature_start-$mature_end";
	if ($mirna ne $prev_mirna) {
		#new precursor miRNA
		push (@mirnas, [$mirna, $accession, $coord, $seq, $mature_coord, $mature_acc, '', '', '', '']);
		$prev_mirna = $mirna;
		$prev_coord = $coord;
	}
	elsif ($prev_coord ne $coord) {
		#for miRNAs expressed from different loci that don't have the -1/-2 name added to the precursor name yet
		$mirnas[-1][2] .= ";$coord";
		$prev_coord = $coord;
	}
	else {
		if ($mature !~ /\*/) {
			#this is the mature strand
			if ($mature_coord ne $mirnas[-1][4]) {
				if ($mirnas[-1][6] eq '') {
					$mirnas[-1][6] = $mature_coord;
					$mirnas[-1][7] = $mature_acc;
				}
				else {
					die "miRNA $mirna has more than 2 mature coords: $mirnas[-1][4], $mirnas[-1][5], $mature_coord\n";
				}
			}
		}
		else {
			#the star strand
			$mirnas[-1][8] = $mature_coord;
			$mirnas[-1][9] = $mature_acc;
		}
	}
}
$sth->finish;
$dbh_mirbase->disconnect;

foreach my $line (@mirnas) {
	print "$line->[0]\t$opt_v\t$line->[1]\t$line->[2]\t$line->[3]\t$line->[4]\t$line->[5]\t$line->[6]\t$line->[7]\t$line->[8]\t$line->[9]\n";
}

sub get_db {
	my $dbname = shift;
	my $dir = dirname(__FILE__);
	my $db_connections = "$dir/../../../config/db_connections.cfg";
	open DB, $db_connections or die "Could not find database connections file $db_connections";
	my @connections = <DB>;
	close DB;
	chomp @connections;
	return split(/\s+/, [grep(/^$dbname/, @connections)]->[0]) or die "Database $dbname not found in $db_connections";
}
