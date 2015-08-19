package Annotate;
require Exporter;
use strict;
use DBI;
use File::Basename;

#common functions for overlapping coordinates in external databases
#used by the scripts that annotate a specific file or annotates all files in a project for clustering

our @ISA = qw(Exporter);
our @EXPORT = qw(build_reference_hashes annotate_coord);

my $BINSIZE = 1000; #size of bin hashes used in coordinate overlap algorithm
#bases past the end of the mature strand that will be considered the stemloop if a tag does not overlap any other explicitly annotated region in the miRNA
my $MIN_OVERLAP = 3; #overlap must be by at least $MIN_OVERLAP bases
my $STEMLOOP_START = 1;
my $STEMLOOP_END = 6;

#reference hashes
my (
	$mirna_co, $mirna_idx,
	$mirna_mature_co, $mirna_mature_idx,
	$mirna_star_co, $mirna_star_idx,
	$mirna_stem_co, $mirna_stem_idx,
	$mirna_pre_co, $mirna_pre_idx,
	$gene_co, $gene_idx,
	$exon_co, $exon_idx,
	$cds_co, $cds_idx,
	$utr5_co, $utr5_idx,
	$utr3_co, $utr3_idx,
	$snorna_co, $snorna_idx,
	$trna_co, $trna_idx,
	$rrna_co, $rrna_idx,
	$snrna_co, $snrna_idx,
	$scrna_co, $scrna_idx,
	$srprna_co, $srprna_idx,
	$rna_co, $rna_idx,
	$rmsk_low_complexity_co, $rmsk_low_complexity_idx,
	$rmsk_ltr_co, $rmsk_ltr_idx,
	$rmsk_line_co, $rmsk_line_idx,
	$rmsk_sine_co, $rmsk_sine_idx,
	$rmsk_satellite_co, $rmsk_satellite_idx,
	$rmsk_simple_repeat_co, $rmsk_simple_repeat_idx,
	$rmsk_dna_co, $rmsk_dna_idx,
	$rmsk_other_co, $rmsk_other_idx,
	$rmsk_unknown_co, $rmsk_unknown_idx
);

sub annotate_coord {
	my ($chr, $start, $end, $strand) = @_;
	$start += $MIN_OVERLAP;
	$end -= $MIN_OVERLAP;
	
	my $overlap_id;
	my %overlap_ids;

	#overlap mature strand, then star strand, the precursor regions, then stemloop, then premature miRNA
	%overlap_ids = overlapcoordinates($mirna_mature_co, $mirna_mature_idx, $chr, $start, $end, $strand);
	if (scalar keys %overlap_ids > 0) {
		($overlap_id) = keys %overlap_ids;
		#mature strand has MIMAT ID in the overlap ID
		my ($mirna, $mimat) = split(',', $overlap_id);
		return ("miRNA", $mirna, "mature,$mimat");
	}
	%overlap_ids = overlapcoordinates($mirna_star_co, $mirna_star_idx, $chr, $start, $end, $strand);
	if (scalar keys %overlap_ids > 0) {
		($overlap_id) = keys %overlap_ids;
		#star strand has MIMAT ID in the overlap ID
		my ($mirna, $mimat) = split(',', $overlap_id);
		return ("miRNA", $mirna, "star,$mimat");
	}
	%overlap_ids = overlapcoordinates($mirna_pre_co, $mirna_pre_idx, $chr, $start, $end, $strand);
	if (scalar keys %overlap_ids > 0) {
		($overlap_id) = keys %overlap_ids;
		return ("miRNA", $overlap_id, "precursor");
	}
	%overlap_ids = overlapcoordinates($mirna_stem_co, $mirna_stem_idx, $chr, $start, $end, $strand);
	if (scalar keys %overlap_ids > 0) {
		($overlap_id) = keys %overlap_ids;
		return ("miRNA", $overlap_id, "stemloop");
	}
	%overlap_ids = overlapcoordinates($mirna_co, $mirna_idx, $chr, $start, $end, $strand);
	if (scalar keys %overlap_ids > 0) {
		($overlap_id) = keys %overlap_ids;
		return ("miRNA", $overlap_id, "unannotated");
	}

	#other small RNAs
	my @smRNA = ('snoRNA', 'tRNA', 'rRNA', 'snRNA', 'scRNA', 'srpRNA', 'rmsk RNA');
	my %smRNA_co = (
		'snoRNA' => $snorna_co,
		'tRNA' => $trna_co,
		'rRNA' => $rrna_co,
		'snRNA' => $snrna_co,
		'scRNA' => $scrna_co,
		'srpRNA' => $srprna_co,
		'rmsk RNA' => $rna_co
	);
	my %smRNA_idx = (
		'snoRNA' => $snorna_idx,
		'tRNA' => $trna_idx,
		'rRNA' => $rrna_idx,
		'snRNA' => $snrna_idx,
		'scRNA' => $scrna_idx,
		'srpRNA' => $srprna_idx,
		'rmsk RNA' => $rna_idx
	);
	foreach my $smRNA (@smRNA) {
		%overlap_ids = overlapcoordinates($smRNA_co{$smRNA}, $smRNA_idx{$smRNA}, $chr, $start, $end, $strand);
		if (scalar keys %overlap_ids > 0) {
			($overlap_id) = keys %overlap_ids;
			return ($smRNA, $overlap_id, "");
		}
	}

	%overlap_ids = overlapcoordinates($gene_co, $gene_idx, $chr, $start, $end, $strand);
	if (scalar keys %overlap_ids > 0) {
		my ($gene) = keys %overlap_ids;
		
		%overlap_ids = overlapcoordinates($exon_co, $exon_idx, $chr, $start, $end, $strand);
		if (scalar keys %overlap_ids > 0) {
			($overlap_id) = keys %overlap_ids;
			return ("Coding Exon", $overlap_id, "No CDS");
		}
		%overlap_ids = overlapcoordinates($utr3_co, $utr3_idx, $chr, $start, $end, $strand);
		if (scalar keys %overlap_ids > 0) {
			($overlap_id) = keys %overlap_ids;
			return ("3 UTR", $overlap_id, "");
		}
		%overlap_ids = overlapcoordinates($utr5_co, $utr5_idx, $chr, $start, $end, $strand);
		if (scalar keys %overlap_ids > 0) {
			($overlap_id) = keys %overlap_ids;
			return ("5 UTR", $overlap_id, "");
		}
		%overlap_ids = overlapcoordinates($cds_co, $cds_idx, $chr, $start, $end, $strand);
		if (scalar keys %overlap_ids > 0) {
			($overlap_id) = keys %overlap_ids;
			return ("Coding Exon", $overlap_id, "");
		}
		return ("Intron", $gene, "");
	}

	#other repeats in UCSC RepeatMasker
	my @repeats = ('LINE', 'SINE', 'LTR', 'Satellite', 'rmsk DNA', 'rmsk Low complexity', 'rmsk Simple repeat', 'rmsk Other', 'rmsk Unknown');
	my %repeats_co = (
		'LINE' => $rmsk_line_co,
		'SINE' => $rmsk_sine_co,
		'LTR' => $rmsk_ltr_co,
		'Satellite' => $rmsk_satellite_co,
		'rmsk DNA' => $rmsk_dna_co,
		'rmsk Low complexity' => $rmsk_low_complexity_co,
		'rmsk Simple repeat' => $rmsk_simple_repeat_co,
		'rmsk Other' => $rmsk_other_co,
		'rmsk Unknown' => $rmsk_unknown_co
	);
	my %repeats_idx = (
		'LINE' => $rmsk_line_idx,
		'SINE' => $rmsk_sine_idx,
		'LTR' => $rmsk_ltr_idx,
		'Satellite' => $rmsk_satellite_idx,
		'rmsk DNA' => $rmsk_dna_idx,
		'rmsk Low idxmplexity' => $rmsk_low_complexity_co,
		'rmsk Simple repeat' => $rmsk_simple_repeat_idx,
		'rmsk Other' => $rmsk_other_idx,
		'rmsk Unknown' => $rmsk_unknown_idx
	);
	foreach my $repeat (@repeats) {
		%overlap_ids = overlapcoordinates($repeats_co{$repeat}, $repeats_idx{$repeat}, $chr, $start, $end, $strand);
		if (scalar keys %overlap_ids > 0) {
			($overlap_id) = keys %overlap_ids;
			return ($repeat, $overlap_id, "");
		}
	}

	return ("Unknown", "", "");
}

sub build_reference_hashes {
	#builds reference hashes of coordinates to do overlaps:
	#mirbase, UCSC knownGene, UCSC rmsk
	#based on algorithm from overlapcoordinates_fast by Richard Corbett

	my $mirbase = shift;
	my $ucsc_db = shift;
	my $species_code = shift;
	my $sample_file = shift; #a sample file in set of files to annotate, to find the labelling formats used

	print STDERR "Building hashes from reference data...\n";

	#first, get the chr formats used in sample file: check if file uses 'chr'chromosome, or just chromosome
	my ($chr_format, $mt_format) = myrna_chr_format($sample_file);

	print STDERR "\tmirBase...";
	($mirna_co, $mirna_idx,
		$mirna_mature_co, $mirna_mature_idx,
		$mirna_star_co, $mirna_star_idx,
		$mirna_stem_co, $mirna_stem_idx,
		$mirna_pre_co, $mirna_pre_idx) = mirbase_hashes($mirbase, $species_code, $chr_format, $mt_format);
	print STDERR "Done\n";

	print STDERR "\tUCSC...\n";
	($gene_co, $gene_idx,
		$exon_co, $exon_idx,
		$cds_co, $cds_idx,
		$utr5_co, $utr5_idx,
		$utr3_co, $utr3_idx,
		$snorna_co, $snorna_idx,
		$trna_co, $trna_idx,
		$rrna_co, $rrna_idx,
		$snrna_co, $snrna_idx,
		$scrna_co, $scrna_idx,
		$srprna_co, $srprna_idx,
		$rna_co, $rna_idx,
		$rmsk_low_complexity_co, $rmsk_low_complexity_idx,
		$rmsk_ltr_co, $rmsk_ltr_idx,
		$rmsk_line_co, $rmsk_line_idx,
		$rmsk_sine_co, $rmsk_sine_idx,
		$rmsk_satellite_co, $rmsk_satellite_idx,
		$rmsk_simple_repeat_co, $rmsk_simple_repeat_idx,
		$rmsk_dna_co, $rmsk_dna_idx,
		$rmsk_other_co, $rmsk_other_idx,
		$rmsk_unknown_co, $rmsk_unknown_idx) = ucsc_hashes($ucsc_db, $chr_format, $mt_format);
	print STDERR "\tDone\nDone\n";
}

sub mirbase_hashes {
	my ($mirbase, $species, $chr_format, $mt_format) = @_;
	my ($dbname, $dbhost, $dbuser, $dbpass) = get_db($mirbase);
	my $dbh_mirbase = DBI->connect("DBI:mysql:database=$dbname;host=$dbhost", $dbuser, $dbpass, {AutoCommit => 0, PrintError => 1}) || die "Could not connect to database: $DBI::errstr";

	#get mirbase species code from organism code
	my $species_code = $dbh_mirbase->selectrow_array("SELECT auto_id FROM mirna_species WHERE organism = '$species'");	
	die "$species organism code not found in database" unless defined $species_code;

	my (%mirna_co, %mirna_idx,
		%mirna_mature_co, %mirna_mature_idx,
		%mirna_star_co, %mirna_star_idx,
		%mirna_stem_co, %mirna_stem_idx,
		%mirna_pre_co, %mirna_pre_idx);

	#complete premature miRNAs
	my ($mirna_id, $chr, $strand, $start, $end);
	my $sth = $dbh_mirbase->prepare("SELECT m.mirna_id, c.xsome, c.contig_start, c.contig_end, c.strand FROM mirna m JOIN mirna_chromosome_build c ON m.auto_mirna = c.auto_mirna WHERE m.auto_species = '$species_code'");
	$sth->execute();
	$sth->bind_columns(\$mirna_id, \$chr, \$start, \$end, \$strand);
	while ($sth->fetchrow_arrayref()) {
		$chr = format_chr($chr, $chr_format, $mt_format);
		addto_overlaphash(\%mirna_co, \%mirna_idx, $chr, $start, $end, $strand, $mirna_id);
	}

	#mature and star strand
	#precursor regions of each
	#estimated stemloop of mature+STEMLOOP_START to mature+STEMLOOP_END
	my ($mature_name, $mature_start, $mature_end, $mature_from, $mature_to, $mature_acc);
	#mature_to and mature_from are in mirna_mature in miRBase v19 and earlier, and moved to mirna_pre_mature in miRBase v20
	$sth = $dbh_mirbase->prepare("SELECT m.mirna_id, ma.mature_name, c.xsome, c.contig_start, c.contig_end, c.strand, mature_from, mature_to, ma.mature_acc
		FROM mirna m 
		JOIN mirna_chromosome_build c ON m.auto_mirna = c.auto_mirna 
		JOIN mirna_pre_mature p ON m.auto_mirna = p.auto_mirna 
		JOIN mirna_mature ma ON p.auto_mature = ma.auto_mature 
		WHERE m.auto_species = '$species_code'");
	$sth->execute();
	$sth->bind_columns(\$mirna_id, \$mature_name, \$chr, \$start, \$end, \$strand, \$mature_from, \$mature_to, \$mature_acc);
	while ($sth->fetchrow_arrayref()) {
		$chr = format_chr($chr, $chr_format, $mt_format);
		if ($strand eq '+') {
			$mature_start = $start + $mature_from - 1;
			$mature_end = $start + $mature_to - 1;
		}
		elsif ($strand eq '-') {
			#"start" is the smaller number, but this corresponds to the 3p end of the miRNA on the - strand
			$mature_start = $end - $mature_to + 1;
			$mature_end = $end - $mature_from + 1;
		}
		else {
			die "Unrecognized miRNA strand: $strand\n";
		}

		#precursor region
		#to determine whether the precursor region is $start to $mature_start or $mature_end to $end, compare ($mature_start - $start) to ($end - $mature_end)
		#1 of those numbers should be over half the miRNA length, while the other is the size of the precursor segment
		my ($pre_start, $pre_end);
		if (($mature_start - $start) < ($end - $mature_end)) {
			$pre_start = $start;
			$pre_end = $mature_start - 1;
		}
		else {
			$pre_end = $end;
			$pre_start = $mature_end + 1;
		}
		addto_overlaphash(\%mirna_pre_co, \%mirna_pre_idx, $chr, $pre_start, $pre_end, $strand, $mirna_id);	

		if ($mature_name =~ /\*/) {
			#star strand
			addto_overlaphash(\%mirna_star_co, \%mirna_star_idx, $chr, $mature_start, $mature_end, $strand, "$mirna_id,$mature_acc");	
		}
		else {
			#mature strand
			addto_overlaphash(\%mirna_mature_co, \%mirna_mature_idx, $chr, $mature_start, $mature_end, $strand, "$mirna_id,$mature_acc");
			my ($stem_start, $stem_end);
			#same logic for determining which arm this strand is as for precursor region
			if (($mature_start - $start) < ($end - $mature_end)) {
				$stem_start = $mature_end + $STEMLOOP_START;
				$stem_end = $mature_end + $STEMLOOP_END;
			}
			else {
				$stem_start = $mature_start - $STEMLOOP_END;
				$stem_end = $mature_start - $STEMLOOP_START;
			}
			addto_overlaphash(\%mirna_stem_co, \%mirna_stem_idx, $chr, $stem_start, $stem_end, $strand, $mirna_id);
		}
	}
	$dbh_mirbase->disconnect;

	return (\%mirna_co, \%mirna_idx, \%mirna_mature_co, \%mirna_mature_idx, \%mirna_star_co, \%mirna_star_idx, \%mirna_stem_co, \%mirna_stem_idx, \%mirna_pre_co, \%mirna_pre_idx);
}

sub ucsc_hashes {
	my ($ucsc_db, $chr_format, $mt_format) = @_;
	my ($dbname, $dbhost, $dbuser, $dbpass) = get_db($ucsc_db);
	my $dbh_ucsc = DBI->connect("DBI:mysql:database=$dbname;host=$dbhost", $dbuser, $dbpass, {AutoCommit => 1, PrintError => 1}) || die "Could not connect to database: $DBI::errstr"; #autocommit set to try to allow auto_reconnect
	$dbh_ucsc->{mysql_auto_reconnect} = 1;
	
	my (%gene_co, %gene_idx,
		%exon_co, %exon_idx,
		%cds_co, %cds_idx,
		%utr5_co, %utr5_idx,
		%utr3_co, %utr3_idx,
		%snorna_co, %snorna_idx,
		%trna_co, %trna_idx,
		%rrna_co, %rrna_idx,
		%snrna_co, %snrna_idx,
		%scrna_co, %scrna_idx,
		%srprna_co, %srprna_idx,
		%rna_co, %rna_idx,
		%rmsk_low_complexity_co, %rmsk_low_complexity_idx,
		%rmsk_ltr_co, %rmsk_ltr_idx,
		%rmsk_line_co, %rmsk_line_idx,
		%rmsk_sine_co, %rmsk_sine_idx,
		%rmsk_satellite_co, %rmsk_satellite_idx,
		%rmsk_simple_repeat_co, %rmsk_simple_repeat_idx,
		%rmsk_dna_co, %rmsk_dna_idx,
		%rmsk_other_co, %rmsk_other_idx,
		%rmsk_unknown_co, %rmsk_unknown_idx);

	#associate repclass name to a pointer to the appropriate hash
	my %repclass_co = (
		'tRNA' => \%trna_co,
		'rRNA' => \%rrna_co,
		'snRNA' => \%snrna_co,
		'scRNA' => \%scrna_co,
		'srpRNA' => \%srprna_co,
		'RNA' => \%rna_co,
		'Low_complexity' => \%rmsk_low_complexity_co,
		'LTR' => \%rmsk_ltr_co,
		'LINE' => \%rmsk_line_co,
		'SINE' => \%rmsk_sine_co,
		'Satellite' => \%rmsk_satellite_co,
		'Simple_repeat' => \%rmsk_simple_repeat_co,
		'DNA' => \%rmsk_dna_co,
		'Other' => \%rmsk_other_co,
		'Unknown' => \%rmsk_unknown_co
	);
	my %repclass_idx = (
		'tRNA' => \%trna_idx,
		'rRNA' => \%rrna_idx,
		'snRNA' => \%snrna_idx,
		'scRNA' => \%scrna_idx,
		'srpRNA' => \%srprna_idx,
		'RNA' => \%rna_idx,
		'Low_complexity' => \%rmsk_low_complexity_idx,
		'LTR' => \%rmsk_ltr_idx,
		'LINE' => \%rmsk_line_idx,
		'SINE' => \%rmsk_sine_idx,
		'Satellite' => \%rmsk_satellite_idx,
		'Simple_repeat' => \%rmsk_simple_repeat_idx,
		'DNA' => \%rmsk_dna_idx,
		'Other' => \%rmsk_other_idx,
		'Unknown' => \%rmsk_unknown_idx
	);

	my ($ref_id, $chr, $strand, $txstart, $txend, $cdsstart, $cdsend, $exonstarts, $exonends);
	my ($exonstarts_arr, $exonends_arr);
	my $sth;

	print STDERR "\t\tsnoRNA...";
	my $ucsc_gene_query = ucsc_gene_table_query($dbh_ucsc, 'snoRNA');
	$sth = $dbh_ucsc->prepare($ucsc_gene_query);
	$sth->execute();
	$sth->bind_columns(\$ref_id, \$chr, \$strand, \$txstart, \$txend, \$cdsstart, \$cdsend, \$exonstarts, \$exonends);
	while ($sth->fetchrow_arrayref()) {
		$chr = format_chr($chr, $chr_format, $mt_format);
		($txstart, $txend, $cdsstart, $cdsend, $exonstarts_arr, $exonends_arr) = format_ucsc($txstart, $txend, $cdsstart, $cdsend, $exonstarts, $exonends);
		my @exonstarts = @{$exonstarts_arr};
		my @exonends = @{$exonends_arr};
		addto_overlaphash(\%gene_co, \%gene_idx, $chr, $txstart, $txend, $strand, $ref_id);
		foreach my $exonindex (0..$#exonstarts) {
			addto_overlaphash(\%snorna_co, \%snorna_idx, $chr, $exonstarts[$exonindex], $exonends[$exonindex], $strand, $ref_id);
		}
	}
	print STDERR "Done\n";

	print STDERR "\t\tKnownGenes with no CDS...";
	my $ucsc_gene_query = ucsc_gene_table_query($dbh_ucsc, 'no_cds');
	$sth = $dbh_ucsc->prepare($ucsc_gene_query);
	$sth->execute();
	$sth->bind_columns(\$ref_id, \$chr, \$strand, \$txstart, \$txend, \$cdsstart, \$cdsend, \$exonstarts, \$exonends);
	while ($sth->fetchrow_arrayref()) {
		$chr = format_chr($chr, $chr_format, $mt_format);
		($txstart, $txend, $cdsstart, $cdsend, $exonstarts_arr, $exonends_arr) = format_ucsc($txstart, $txend, $cdsstart, $cdsend, $exonstarts, $exonends);
		my @exonstarts = @{$exonstarts_arr};
		my @exonends = @{$exonends_arr};
		addto_overlaphash(\%gene_co, \%gene_idx, $chr, $txstart, $txend, $strand, $ref_id);
		foreach my $exonindex (0..$#exonstarts) {
			addto_overlaphash(\%exon_co, \%exon_idx, $chr, $exonstarts[$exonindex], $exonends[$exonindex], $strand, $ref_id);
		}
	}
	print STDERR "Done\n";

	print STDERR "\t\tKnownGenes...";
	$ucsc_gene_query = ucsc_gene_table_query($dbh_ucsc, 'cds');
	$sth = $dbh_ucsc->prepare($ucsc_gene_query);
	$sth->execute();
	$sth->bind_columns(\$ref_id, \$chr, \$strand, \$txstart, \$txend, \$cdsstart, \$cdsend, \$exonstarts, \$exonends);
	while ($sth->fetchrow_arrayref()) {
		$chr = format_chr($chr, $chr_format, $mt_format);
		($txstart, $txend, $cdsstart, $cdsend, $exonstarts_arr, $exonends_arr) = format_ucsc($txstart, $txend, $cdsstart, $cdsend, $exonstarts, $exonends);
		my @exonstarts = @{$exonstarts_arr};
		my @exonends = @{$exonends_arr};
		addto_overlaphash(\%gene_co, \%gene_idx, $chr, $txstart, $txend, $strand, $ref_id);

		#split exons into coding, 3'utr and 5'utr
		my ($utr5_start, $utr5_end, $utr3_start, $utr3_end);
		if ($strand eq '+') {
			$utr5_start = $exonstarts[0];
			$utr5_end = $cdsstart - 1;
			$utr3_start = $cdsend + 1;
			$utr3_end = $exonends[-1];
		}
		elsif ($strand eq '-') {
			#on the minus strand, 5 and 3 are reversed
			$utr5_start = $cdsend + 1;
			$utr5_end = $exonends[-1];
			$utr3_start = $exonstarts[0];
			$utr3_end = $cdsstart - 1;
		}
		else {
			die "Strand $strand not recognized for UCSC gene $ref_id";
		}
		addto_overlaphash(\%utr5_co, \%utr5_idx, $chr, $utr5_start, $utr5_end, $strand, $ref_id);
		addto_overlaphash(\%utr3_co, \%utr3_idx, $chr, $utr3_start, $utr3_end, $strand, $ref_id);
		foreach my $exonindex (0..$#exonstarts) {
			my $exonstart = $exonstarts[$exonindex];
			my $exonend = $exonends[$exonindex];
			next if $exonend < $cdsstart || $exonstart > $cdsend; #skip exon if it's completely outside CDS (shouldn't happen)
			$exonstart = $cdsstart if $exonstart < $cdsstart;
			$exonend = $cdsend if $exonend > $cdsend;
			addto_overlaphash(\%cds_co, \%cds_idx, $chr, $exonstart, $exonend, $strand, $ref_id);
		}
	}
	print STDERR "Done\n";

	print STDERR "\t\tRepeatMasker...";
	my ($start, $end, $repclass);
	my @rmsk_chr_tables = ucsc_rmsk_chr_tables($dbh_ucsc);
	my %new_classes; #warn about new repeat classes not tracked by scripts
	foreach my $chr_table (@rmsk_chr_tables) {
		$sth = $dbh_ucsc->prepare("SELECT genoName, genoStart, genoEnd, strand, repName, repClass	FROM $chr_table");
		$sth->execute();
		$sth->bind_columns(\$chr, \$start, \$end, \$strand, \$ref_id, \$repclass);
		while ($sth->fetchrow_arrayref()) {
			unless (exists $repclass_co{$repclass}) {
				$new_classes{$repclass}++;
				next;
			}
			$chr = format_chr($chr, $chr_format, $mt_format);
			addto_overlaphash($repclass_co{$repclass}, $repclass_idx{$repclass}, $chr, $start, $end, $strand, $ref_id);
		}
	}
	foreach my $new_class (keys %new_classes) {
		warn "Warning: New repeat class $new_class found, $new_classes{$new_class} instances in database.\n";
	}
	print STDERR "Done\n";
	$dbh_ucsc->disconnect;
	
	return (\%gene_co, \%gene_idx,
		\%exon_co, \%exon_idx,
		\%cds_co, \%cds_idx,		
		\%utr5_co, \%utr5_idx,
		\%utr3_co, \%utr3_idx,
		\%snorna_co, \%snorna_idx,
		\%trna_co, \%trna_idx,
		\%rrna_co, \%rrna_idx,
		\%snrna_co, \%snrna_idx,
		\%scrna_co, \%scrna_idx,
		\%srprna_co, \%srprna_idx,
		\%rna_co, \%rna_idx,
		\%rmsk_low_complexity_co, \%rmsk_low_complexity_idx,
		\%rmsk_ltr_co, \%rmsk_ltr_idx,
		\%rmsk_line_co, \%rmsk_line_idx,
		\%rmsk_sine_co, \%rmsk_sine_idx,
		\%rmsk_satellite_co, \%rmsk_satellite_idx,
		\%rmsk_simple_repeat_co, \%rmsk_simple_repeat_idx,
		\%rmsk_dna_co, \%rmsk_dna_idx,
		\%rmsk_other_co, \%rmsk_other_idx,
		\%rmsk_unknown_co, \%rmsk_unknown_idx);
}

sub ucsc_gene_table_query {
	#used in build_reference_hashes
	#For the mirror of the genes table, check if knownGenes is available
	#if not, use refGenes
	#if neither of those are available, die
	my $dbh_ucsc = shift;
	my $gene_type = shift; #one of no_cds, cds, snoRNA

	my $subclause = "";
	if ($gene_type eq 'snoRNA') {
		$subclause = " WHERE cdsStart = cdsEnd AND (x.geneSymbol LIKE 'SNOR%' OR x.description LIKE '% small nucleolar RNA%')";
	}
	elsif ($gene_type eq 'no_cds') {
		$subclause = " WHERE cdsStart = cdsEnd AND NOT (x.geneSymbol LIKE 'SNOR%' OR x.description LIKE '% small nucleolar RNA%')";
	}
	else {
		$subclause = " WHERE cdsStart != cdsEnd";
	}

	my $ucsc_gene_table;
	my $ucsc_gene_query;	
	my $sth = $dbh_ucsc->prepare('SHOW TABLES LIKE ?');
	$sth->execute('knownGene');
	$sth->bind_columns(\$ucsc_gene_table);
	$sth->fetch;
	$ucsc_gene_query = "SELECT x.geneSymbol AS name, kg.chrom, kg.strand, kg.txStart, kg.txEnd, kg.cdsStart, kg.cdsEnd, kg.exonStarts, kg.exonEnds FROM knownGene kg JOIN kgXref x ON x.kgID = kg.name".$subclause;
	unless (defined $ucsc_gene_table) {
		$sth->execute('refGene');
		$sth->bind_columns(\$ucsc_gene_table);
		$sth->fetch;
		$ucsc_gene_query = "SELECT x.geneSymbol AS name, r.chrom, r.strand, r.txStart, r.txEnd, r.cdsStart, r.cdsEnd, r.exonStarts, r.exonEnds FROM refGene r JOIN kgXref x ON x.refseq = r.name".$subclause;
	}
	die "Neither knownGene nor refGene table found in UCSC database, update script to find another source for genes\n" unless defined $ucsc_gene_table;
	$sth->finish;
	return $ucsc_gene_query;
}

sub ucsc_rmsk_chr_tables {
	#used in build_reference_hashes
	#checks if the rmsk table for the organism is 1 table or divided by chr
	my $dbh_ucsc = shift;
	my (@chr_tables, $chr_table);
	my $sth = $dbh_ucsc->prepare('SHOW TABLES LIKE ?');
	$sth->execute('%rmsk');
	$sth->bind_columns(\$chr_table);
	while ($sth->fetchrow_arrayref()) {
		#if there is more than 1 underscore (to separate chr from rmsk), don't count the table
		my $num_ = ($chr_table =~ tr/_/_/);
		next if $num_ > 1;
		push(@chr_tables, $chr_table);
	}
	$sth->finish;
	return @chr_tables;
}

sub addto_overlaphash {
	#used in build_reference_hashes
	#takes a hash for overlap coordinates, by reference, and adds the coordinate sets to the appropriate places
	my $loc_ref = shift;
	my $list_ref = shift;
	my ($chr, $start, $end, $strand, $id) = @_;
	my $chrstrand = $chr.$strand;

	while (exists($loc_ref->{$id})) {
		#if $id already exists for a set of coordinates, append ___num after it to make a unique key for the hash
		#this method keeps trying to append ___nums until it hits a unique one
		my ($num) = $id =~ /___(\d+)/;
		$id = $id.'___1' unless defined $num;
		$num = int(rand(1000000000000000));
		$id =~ s/___\d+/___$num/;
	}

	#save the info for this element in a hash
  $loc_ref->{$id}{set}{"$start-$end"} = 1;
  $loc_ref->{$id}{chr} = $chrstrand;

  #Add this element to a hash table to speed up the lookups later on
  my $start_ind = int($start/$BINSIZE);
  my $stop_ind = int($end/$BINSIZE);
  foreach( $start_ind-1 ..$stop_ind+1 ) {
    $list_ref->{$chrstrand}{$_}{$id} = 1;
  }
}

sub overlapcoordinates {
	my $loc_ref = shift;
	my $list_ref = shift;
	my ($chr, $start, $end, $strand, $olap_all) = @_;
	my $chrstrand = $chr.$strand;
	$olap_all = 0 unless defined $olap_all;
	#if olap_all is true, overlap goes through and finds every match, otherwise return the first match found
	my %return_hash;

  #find all the features that overlap with the given region
  my $strt = int($start/$BINSIZE);
  my $stp = int($end/$BINSIZE);
  foreach( $strt ..$stp ) {
		next unless exists $list_ref->{$chrstrand} && exists $list_ref->{$chrstrand}{$_};
    foreach my $ref_id (keys %{ $list_ref->{$chrstrand}{$_} } ) {
      #if there is any overlap
			foreach my $coords (keys %{$loc_ref->{$ref_id}{set}}) {
				my ($ref_start, $ref_end) = split('-', $coords);
				if ($ref_start <= $end && $ref_end >= $start) {
					#remove ___\d+ added when building the hash
					$ref_id =~ s/___\d+//;
					$return_hash{$ref_id} = 1	;			
					return %return_hash unless $olap_all;
				}
			}
    }
  }
	return %return_hash;
}

sub myrna_chr_format {
	#used in build_reference_hashes
	#checks chromosome naming format in myRNA so that coordinates retrieved in reference databases can be formatted the same way to make comparisons possible
	my $sample_file = shift;
	print STDERR "\tChecking myRNA formatting to format reference data...";
	#errs on the side of not having 'chr'
	#for .sam files, look for the @SQ header and sees if there is a chr in there
	#for custom files without headers or unaligned regions (peaks from FindPeaks, pulls out the first line and looks for 'chr'
	my $has_chr = `grep -m1 "^\@SQ" $sample_file | grep -c chr`;
	$has_chr = $has_chr || `head -1 $sample_file | grep -c chr`;
	chomp $has_chr;
	my $chr_format = $has_chr ? 'chr' : '';
	my $has_mt = `grep -c -m1 ${chr_format}MT $sample_file`;
	chomp $has_mt;
	my $mt_format = $has_mt ? 'MT' : 'M';
	print STDERR "Done\n";
	return $chr_format, $mt_format;
}

sub format_chr {
	#used in build_reference_hashes
	#formats given chr to supplied format
	my $chr = shift;
	my $chr_format = shift;
	my $mt_format = shift;

	#if myRNA database uses 'chr1' format, format has 'chr'
	if ($chr_format eq 'chr' && $chr !~ /chr/) {
		$chr = 'chr'.$chr;
	}
	#if myRNA does not use 'chr' format, or it's chromosomes don't have 'chr' (eg. scaffold), format doesn't have 'chr'
	elsif ($chr_format eq '') {
		$chr =~ s/^chr//;
	}

	#change chrM or MT to match myRNA
	if ($chr eq $chr_format.'M' || $chr eq $chr_format.'MT') {
		$chr = $chr_format.$mt_format;
	}
	return $chr;
}

sub format_ucsc {
	my ($txstart, $txend, $cdsstart, $cdsend, $exonstarts, $exonends) = @_;
	my @exonstarts = split(',', $exonstarts);
	my @exonends = split(',', $exonends);

	#add 1 to 0 based ucsc coordinates
	$txstart += 1;
	$txend += 1;
	$cdsstart += 1;
	$cdsend += 1;
	foreach my $exonindex (0..$#exonstarts) {
		$exonstarts[$exonindex] += 1;
		$exonends[$exonindex] += 1;
	}

	return ($txstart, $txend, $cdsstart, $cdsend, \@exonstarts, \@exonends);
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

1;
