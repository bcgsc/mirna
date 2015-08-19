#!/usr/bin/env perl
use strict; 
use Getopt::Std; 
use vars qw ($opt_m $opt_o $opt_p $opt_h);
getopts('m:o:p:h:');
use DBI;
use File::Find;
use File::Basename;

my $help_text = "
############## $0 ########################

Using a list of genes and a set of expression files, creates an expression matrix.
Uses miRNA.txt and crossmapped.txt files created by alignment_stats.pl.

3 versions of the expression matrix are created.
expn_matrix_mimat.txt: matrix with raw read counts
expn_matrix_mimat_norm.txt: matrix with normalized counts (counts per million mature miRNA aligned tags)
expn_matrix_mimat_norm_log.txt: as above, but values are taken to log2

Usage: $0 -m mirbase_db -o species_code -p project_dir

########################################################
";

if ($opt_h || !($opt_m && $opt_o && $opt_p)) { print $help_text; exit 1; }

#Data structure
my %genes; #hash gene->sample->value
my @genes; #full list of gene names
my %samples; #hash of samples->total tags
my %sample_names; #sample names
my %mimats; #maps mimat id to a unique mirna id

my $norm_factor = 1000000; #Normalize to this number of tags

print STDERR "Searching for input files in project directory...\n";
my @mirnafiles;
find(\&findinputfiles, $opt_p);
print STDERR "Done\n";
die "No input files found in $opt_p" unless scalar @mirnafiles;
@mirnafiles = sort @mirnafiles;

#Read in gene list
my ($genes, $mimats) = get_mirna_genes($opt_m, $opt_o);
@genes = @$genes;
%mimats = %$mimats;

#Read in mirna expression
foreach my $source (@mirnafiles) {
	my $sample = (split(/\//, dirname($source)))[-1];
	$sample =~ s/_features//;
	$sample_names{$sample} = 1;
	
	my $sample_total = 0; #Total tag count for sample based on expr file
	open(FH, $source) || die "Cannot open expression file $source: $!";
	while(<FH>) {
		chomp;
		my ($mirna, $value) = split(' ');
		my @mir_entries = split(';', $mirna); #crossmapped miRNAs are split by ';', if the miRNA is not crossmapped, mir_entries will only be 1 element equal to the original $mirna
		foreach my $mir_entry (@mir_entries) {
			my ($gene, $mature, $acc) = split(',', $mir_entry);
			next unless $mature eq 'mature' || $mature eq 'star'; #skip precursor, stemloop, and unannotated
			my $id = $mimats{$acc}.".".$acc;
			$genes{$id}{$sample} += $value;
			$sample_total += $value;
		}
	}
	close(FH);
	$samples{$sample} += $sample_total;
}

#Header
my $header = "Gene";
foreach my $sample (sort keys %sample_names) {
	$header .= "\t$sample";
}
$header .= "\n";

#Data
open RAW, ">$opt_p/expn_matrix_mimat.txt" or die "Can't write out expression matrix $opt_p/expn_matrix_mimat.txt: $!\n";
open NORM, ">$opt_p/expn_matrix_mimat_norm.txt" or die "Can't write out expression matrix $opt_p/expn_matrix_mimat_norm.txt: $!\n";
open LOG, ">$opt_p/expn_matrix_mimat_norm_log.txt" or die "Can't write out expression matrix $opt_p/expn_matrix_mimat_norm_log.txt: $!\n";
print RAW $header;
print NORM $header;
print LOG $header;

foreach my $gene (@genes) {
	print RAW "$gene";
	print NORM "$gene";
	print LOG "$gene";

	foreach my $sample (sort keys %sample_names) {
		my $count = $genes{$gene}{$sample} || 0;
		warn "No reads mapped to miRNAs found in sample $sample. Perhaps the mirna_species.txt file is empty?\n" if $samples{$sample} == 0;
		my $norm_count;
		my $log_count;
		if ($samples{$sample} == 0) {
			$norm_count = 0;
			$log_count = 0;
		}
		else {
			$norm_count = sprintf("%.6f", $count*$norm_factor/$samples{$sample});
			$log_count = sprintf("%.6f", $norm_count == 0 ? 0 : log($norm_count) / log(2));
		}

		print RAW "\t$count";
		print NORM "\t$norm_count";
		print LOG "\t$log_count";
	}
	print RAW "\n";
	print NORM "\n";
	print LOG "\n";
}

close RAW;
close NORM;
close LOG;
	
sub get_mirna_genes {
	my $db = shift;
	my $species = shift;
	my ($dbname, $dbhost, $dbuser, $dbpass) = get_db($db);
	my $dbh_mirbase = DBI->connect("DBI:mysql:database=$dbname;host=$dbhost", $dbuser, $dbpass, {AutoCommit => 0, PrintError => 1}) || die "Could not connect to database: $DBI::errstr";
	#get mirbase species code from organism code
	my $species_code = $dbh_mirbase->selectrow_array("SELECT auto_id FROM mirna_species WHERE organism = '$species'");
	die "$species_code organism code not found in database" unless defined $species_code;
	
	my (%mirnas, %mimats);
	my ($mirna_id, $mature_acc);
	#only use miRNAs where there is coordinate information and mature strand information
	my $sth = $dbh_mirbase->prepare("SELECT DISTINCT mirna_id, ma.mature_acc FROM mirna m JOIN mirna_chromosome_build co ON m.auto_mirna = co.auto_mirna JOIN mirna_pre_mature p ON m.auto_mirna = p.auto_mirna JOIN mirna_mature ma ON p.auto_mature = ma.auto_mature WHERE auto_species = '$species_code' ORDER BY mirna_id");	
	$sth->execute();
	$sth->bind_columns(\$mirna_id, \$mature_acc);
	while ($sth->fetchrow_arrayref()) {
		#if this MIMAT has not been seen yet, just add the mirna.mimat to the hash
		unless (exists $mimats{$mature_acc}) {
			$mimats{$mature_acc} = $mirna_id; #track the miRNA ID associated with the MIMAT (explained in trim_id)
			$mirnas{"$mirna_id.$mature_acc"} = 1;
			next;
		}
		my $trimmed_mirna_id = trim_id($mirna_id, $mimats{$mature_acc}); #find common miRNA_id to use for MIMAT
		delete($mirnas{"$mimats{$mature_acc}.$mature_acc"});
		$mirnas{"$trimmed_mirna_id.$mature_acc"} = 1;
		$mimats{$mature_acc} = $trimmed_mirna_id;
	}
	$sth->finish;
	$dbh_mirbase->disconnect;
	my @mirnas = sort keys %mirnas;
	return \@mirnas, \%mimats;
}

sub trim_id {
	#if this MIMAT has been seen before, it belongs to multiple similar miRNAs
	#use a "common" miRNA name which is the common component shared between miRNA names
	#eg. 
	#hsa-mir-181a and hsa-mir-181b gets trimmed down to hsa-mir-181
	#hsa-mir-941-1 and hsa-mir-941-2 gets trimmed down to hsa-mir-941
	#hsa-mir-941 and hsa-mir-941-3 gets trimmed down to hsa-mir-941 (this happens when the miRNA_id has been processed once already)
	my $id1 = shift;
	my $id2 = shift;

	my $minlen = length($id1);
	$minlen = length($id2) if length($id2) < $minlen;

	#walk backwards through the strings until they're equal
	for (my $i = $minlen; $i >=0; $i--) {
		my $sub_id1 = substr($id1, 0, $i);
		my $sub_id2 = substr($id2, 0, $i);
		if ($sub_id1 eq $sub_id2) {
			#remove trailing - if necessary
			$sub_id1 =~ s/-$//;
			return $sub_id1;
		}
	}
	die "No common substring between $id1 and $id2 for the same MIMAT."
}

sub findinputfiles {
	if (($File::Find::name =~ /_features\/miRNA.txt$/ || $File::Find::name =~ /_features\/crossmapped.txt$/ ) && $File::Find::name !~ /obsoleted/) {
		#do not look in the standard "obsoleted" subfolder
		push(@mirnafiles, $File::Find::name);
		print STDERR "\t$File::Find::name\n";
	}
}

sub get_db {
	my $dbname = shift;
	my $dir = dirname(__FILE__);
	my $db_connections = "$dir/../../config/db_connections.cfg";
	open DB, $db_connections or die "Could not find database connections file $db_connections";
	my @connections = <DB>;
	close DB;
	chomp @connections;
	return split(/\s+/, [grep(/^$dbname/, @connections)]->[0]) or die "Database $dbname not found in $db_connections";
}
