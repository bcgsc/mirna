#!/usr/bin/env perl
use strict; 
use Getopt::Std; 
use vars qw ($opt_m $opt_o $opt_p $opt_h);
getopts('m:o:p:h:');
use File::Find;
use File::Basename;

my $help_text = "
############## $0 ########################

Using a list of genes and a set of expression files, creates an expression matrix.

3 versions of the expression matrix are created.
expn_matrix_mimat.txt: matrix with raw read counts
expn_matrix_mimat_norm.txt: matrix with normalized counts (counts per million mature miRNA aligned tags)
expn_matrix_mimat_norm_log.txt: as above, but values are taken to log2 (log(0) is written out as 0)

Usage: $0 -m miRNA_ADF -p Level_3_archive_directory 

########################################################
";

if ($opt_h || !($opt_m && $opt_p)) { print $help_text; exit 1; }

my $norm_factor = 1000000; #Normalize to this number of tags
my $log0 = "0"; #log(0) can be written out as other strings, eg. "NA"

#Data structure
my %genes; #hash gene->sample->value
my @genes; #full list of gene names
my %samples; #hash of samples->total tags
my %sample_names; #sample names
my %mimats; #maps mimat id to a unique mirna id

print STDERR "Searching for input files in project directory...\n";
my @mirnafiles;
find(\&findinputfiles, $opt_p);
print STDERR "Input files found\n";
die "No input files found in $opt_p" unless scalar @mirnafiles;
@mirnafiles = sort @mirnafiles;

#Read in gene list
my ($genes, $mimats) = get_mirna_genes($opt_m);
@genes = @$genes;
%mimats = %$mimats;

#Read in mirna expression
foreach my $source (@mirnafiles) {
	my $sample = basename($source, ".isoform.quantification.txt");
	$sample_names{$sample} = 1;
	
	my $sample_total = 0; #Total tag count for sample based on expr file
	open(FH, $source) || die "Cannot open expression file $source: $!";
	while(<FH>) {
		chomp;
		my ($mirna, $coords, $count, $norm_count, $crossmap, $annot) = split;
		#ignore anything that's not a mature or star strand
		#count cross-mapped reads towards each annotation they map to
		my ($mature, $acc) = split(',', $annot);
		next unless $mature eq 'mature' || $mature eq 'star';
		my $id = $mimats{$acc}.".".$acc;
		$genes{$id}{$sample} += $count;
		$sample_total += $count;
	}
	close(FH);
	$samples{$sample} = $sample_total;
}

#Expression matrix header
my $header = "Gene";
foreach my $sample (sort keys %sample_names) {
	$header .= "\t$sample";
}
$header .= "\n";

#Write expression matrix
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
		my $norm_count;
		my $log_count;
		if ($samples{$sample} == 0) {
			$norm_count = 0;
			$log_count = $log0;
		}
		else {
			$norm_count = sprintf("%.6f", $count*$norm_factor/$samples{$sample});
			$log_count = sprintf("%.6f", $norm_count == $log0 ? 0 : log($norm_count) / log(2));
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
	my $adf = shift;
	
	my (%mirnas, %mimats);

	open FH, "<$adf" or die "Cannot read ADF reference $adf: $!\n";
	while (<FH>) {
		my @fields = split;
		chomp @fields;
		#adf format:
		#miRNA_ID miRBase_Ver Accession Genomic_Coords Precursor_Seq Mature_Coords Mature_Accession Alt_Mature_Coords Alt_Mature_Accession Star_Coords Star_Accession
		#sorted by miRNA_ID
		my $mirna_id = $fields[0];
		next if $mirna_id eq "miRNA_ID"; #skip header
		my @mimats = ($fields[6]);
		push(@mimats, $fields[8]) unless $fields[8] eq "";
		push(@mimats, $fields[10]) unless $fields[10] eq "";
		foreach my $mimat (@mimats) {
			#if this MIMAT has not been seen yet, just add the mirna.mimat to the hash
			unless (exists $mimats{$mimat}) {
				$mimats{$mimat} = $mirna_id; #track the miRNA ID associated with the MIMAT (explained in trim_id)
				$mirnas{"$mirna_id.$mimat"} = 1;
				next;
			}
			my $trimmed_mirna_id = trim_id($mirna_id, $mimats{$mimat}); #find common miRNA_id to use for MIMAT
			delete($mirnas{"$mimats{$mimat}.$mimat"});
			$mirnas{"$trimmed_mirna_id.$mimat"} = 1;
			$mimats{$mimat} = $trimmed_mirna_id;
		}
	}
	close FH;
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
	if ($File::Find::name =~ /.isoform.quantification.txt$/) {
		push(@mirnafiles, $File::Find::name);
		print STDERR "\t$File::Find::name\n";
	}
}
