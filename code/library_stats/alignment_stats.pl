#!/usr/bin/env perl
use strict;
use Getopt::Std;
use Pod::Usage;
use File::Find;
use File::Basename;

use vars qw($opt_p);
getopts("p:");

#finds all .sam files in a given base directory, assuming they've been annotated with XC
#also assumes sam file names are $library[_$index].sam
#produce alignment stats for each library
#produce a set of strand+.bed and strand-.bed files for peak discovery used in novel miRNA prediction

my $usage = "$0 -p project_directory\n";
die "$usage" unless $opt_p;

print STDERR "Searching for .sam files in project directory...\n";
my @samfiles;
find(\&findsamfiles, $opt_p);
print STDERR "Done\n";
die "No sam files found in $opt_p" unless scalar @samfiles;
@samfiles = sort @samfiles;

my $MAX_XA = 10; #number of alignments a read can have, as run by `bwa samse -n` by LIMS
my $MAX_NM = 0; #most mismatches allowed in the alignment, according to NM flag
my $MAX_X0 = 3; #most alignments a read can have before it's ignored
#annotation priority from XA and XD terms in sam file
#uses XD terms where available, otherwise XC terms
my %ANNOTATION_PRIORITY = (
	#mirnas are priorities 0 to 9
	'mature' => 0,
	'star' => 1,
	'precursor' => 2,
	'stemloop' => 3,
	'unannotated' => 4,

	#priorities 10 to 90 are other possible annotations
	#snoRNAs from UCSC genes
	'snoRNA' => 11,

	#small RNAs from UCSC RepeatMasker
	'tRNA' => 15,
	'rRNA' => 16,
	'snRNA' => 17,
	'scRNA' => 18,
	'srpRNA' => 19,
	'rmsk RNA' => 20,
	
	#genes
	'No CDS' => 30, #genes with no CDS regions in UCSC are likely various types of RNAs
	'3 UTR' => 31,
	'5 UTR' => 32,
	'Coding Exon' => 33,
	'Intron' => 34,

	#other elements in UCSC RepeatMasker
	'LINE' => 50,
	'SINE' => 51,
	'LTR' => 52,
	'Satellite' => 53,
	'rmsk DNA' => 54,
	'rmsk Low complexity' => 55,
	'rmsk Simple repeat' => 56,
	'rmsk Other' => 57,
	'rmsk Unknown' => 58,

	#priorities > 90 have no species details
	'Unknown' => 99,
	#priorities >= 100 aren't included in any output
	'' => 100
);
#list of features in "human readable" format, for outputting to CSV
my @feature_list = ("mature miRNA","* miRNA","precursor miRNA","miRNA loop","unannotated miRNA","snoRNA","tRNA", 'rRNA', 'snRNA', 'scRNA', 'srpRNA', 'Other RepeatMasker RNAs', "RNA (No CDS)","3' UTR","5' UTR","Coding Exon","Intron","LINE","SINE","LTR","Satellite","RepeatMasker DNA", "RepeatMasker Low complexity", "RepeatMasker Simple repeat", "RepeatMasker Other", "RepeatMasker Unknown", "Unknown");
#keys in %ANNOTATION_PRIORITY, sorted in order of @feature_list
my @feature_key_list = ('mature', 'star', 'precursor', 'stemloop', 'unannotated', 'snoRNA', 'tRNA', 'rRNA', 'snRNA', 'scRNA', 'srpRNA', 'rmsk RNA', 'No CDS', '3 UTR', '5 UTR', 'Coding Exon', 'Intron', 'LINE', 'SINE', 'LTR', 'Satellite', 'rmsk DNA', 'rmsk Low complexity', 'rmsk Simple repeat', 'rmsk Other', 'rmsk Unknown', 'Unknown');

#adapter positions to group by for reporting adapter trimming summary
my @ADAPTER_CATEGORIES = (0, 14, 25, 35, -1); #-1 length represents infinite length
#labels for csv
my @adapter_names = ("Adapter dimers", "Adapter at 1-14bp", "Adapter at 15-25bp", "Adapter at 26-35bp", "Adapter after 35bp");

my %BED; #hashes for file handles to write to miRNA, not_miRNA, and unknown bed files for peak discovery, in the format $BED{category}{chr}{strand}
my $beddir;

open PROJ, ">$opt_p/alignment_stats.csv" or die "Can't write out alignment stats file: $opt_p/alignment_stats.csv: $!\n";
print PROJ "Library,Index,Total Reads >= 15bp,";
print PROJ "% Adapter dimers,";
print PROJ join(',', @adapter_names);
print PROJ ",Aligned Reads Post-Filter,% Aligned Reads,Unaligned Reads,% Unaligned Reads,Filtered Reads without XA,Softclipped Reads,Chastity Failed Reads Post-Filter,miRNA Species,miRNA Species Covered by >= 10 Reads";
print PROJ ",Total miRNA,Crossmapped miRNA,".join(',', @feature_list);
print PROJ ",% Total miRNA,% Crossmapped miRNA,% ".join(',% ', @feature_list);
print PROJ "\n";

foreach my $samfile (@samfiles) {
	my $samdir = dirname($samfile);
	my $filename = basename($samfile);
	$filename =~ s/\.sam$//;
	my ($lib, $index) = split('_', $filename);
	$index = '' unless defined $index;

	my @adapter_counts = read_adapter($samdir."/".$filename."_adapter.report");
	my $total_reads = 0; #counts up all reads in sample, including those not sent for alignment due to adapter
	foreach my $adapter_count (@adapter_counts) {
		$total_reads += $adapter_count;
	}

	my $outdir = "$samdir/$filename\_features";
	unless (-d $outdir) {
		mkdir $outdir or die "Could not create subdirectory for feature files $outdir: $!";
	}
	$beddir = "$outdir/bed";
	unless (-d $beddir) {
		mkdir $beddir or die "Could not create subdirectory for bed files $beddir: $!";
	}

	my ($feature_species, $mirna_species, $filtered_taglen, $crossmapped, $num_reads, $num_filtered, $num_unaligned, $num_overfiltered, $num_softclipped, $num_chastity, $softclip_features, $chastity_taglen, $isoforms) = process_file($samfile);

	#add line to csv
	my $mirna_diversity = scalar keys %{$mirna_species};
	my $pct_aligned = sprintf('%.2f%%', $num_filtered * 100 / $num_reads);
	my $pct_unaligned = sprintf('%.2f%%', $num_unaligned * 100 /$num_reads);
	my $pct_adapter_dimer = sprintf('%.2f%%', $total_reads == 0 ? 0 : $adapter_counts[0] * 100 /$total_reads);
	print PROJ "$lib,$index,$num_reads,$pct_adapter_dimer,".join(',', @adapter_counts).",$num_filtered,$pct_aligned,$num_unaligned,$pct_unaligned,$num_overfiltered,$num_softclipped,$num_chastity,$mirna_diversity,";

	#create a hash of miRNAs without annotation, to find read coverage for each miRNA locus
	my %mirna_locus;
	my $mirna_diversity_filtered = 0; #counts of miRNAs supported by at least 10 reads anywhere in the entire locus
	foreach my $mirna_annotation (keys %{$feature_species->{miRNA}}) {
		my $key = (split(',', $mirna_annotation))[0];
		$mirna_locus{$key} += $feature_species->{miRNA}{$mirna_annotation};
	}
	#add crossmapped reads into mirna_locus as well, the reads that are crossmapped (as opposed to the uniquely mapped reads that map to a crossmapped miRNA), will be counted multiple times
	#at the same time, count up the number of reads in crossmapped miRNAs total
	my $num_crossmapped = 0;
	foreach my $crosslabel (keys %{$crossmapped}) {
		$num_crossmapped += $crossmapped->{$crosslabel};

		my @mirna_annotations = split(';', $crosslabel); #break up the multiple miRNAs in the label
		foreach my $mirna_annotation (@mirna_annotations) {
			my $mirna = (split(',', $mirna_annotation))[0];
			$mirna_locus{$mirna} += $crossmapped->{$crosslabel};
		}
	}
	
	foreach my $mirna (keys %mirna_locus) {
		$mirna_diversity_filtered++ if $mirna_locus{$mirna} >= 10;
	}
	print PROJ "$mirna_diversity_filtered,";

	#for feature counts, run through feature_species to get features that have not been crossmapped
	my (%num_features);
	my $num_mirna = 0; #sum up all miRNA reads in crossmapped + ANNOTATION_PRIORITY < 10
	foreach my $xc (keys %{$feature_species}) {
		foreach my $species (keys %{$feature_species->{$xc}}) {
			my ($xi, $xd) = split(',', $species);
			my $annot_key = get_annot_key($xc, $xd);
			$num_features{$annot_key} += $feature_species->{$xc}{$species};
			$num_mirna += $feature_species->{$xc}{$species} if $ANNOTATION_PRIORITY{$annot_key} < 10;
		}
	}

	#for features that don't have individual species, use filtered_taglen to get the read counts
	for my $feature (@feature_key_list) {
		if ($ANNOTATION_PRIORITY{$feature} > 90) {
			foreach my $taglength (keys %{$filtered_taglen}) {
				$num_features{$feature} += $filtered_taglen->{$taglength}{$feature};
			}
		}
	}

	$num_mirna += $num_crossmapped;
	print PROJ "$num_mirna,$num_crossmapped,";
	foreach my $feature (@feature_key_list) {
		$num_features{$feature} = 0 unless exists $num_features{$feature};
		print PROJ "$num_features{$feature},";
	}
	printf PROJ "%.2f%%,", $num_mirna * 100 / $num_filtered;
	printf PROJ "%.2f%%,", $num_crossmapped * 100 / $num_filtered;
	my $str; #sprintf the results to a string, then strip off the final comma before printing to file
	foreach my $feature (@feature_key_list) {
		$str .= sprintf("%.2f%%,", $num_features{$feature} * 100 /$num_filtered);
	}
	chop $str;
	print PROJ "$str\n";
	
	#write out feature files listing species coverage
	write_feature_species($feature_species, $outdir, $num_mirna);
	#write out crossmapped miRNAs
	write_crossmapped($crossmapped, $outdir, $num_mirna);
	#write out miRNA locus read counts
	write_mirna_locus(\%mirna_locus, $outdir, $num_mirna);
	#write out filtered tag lengths
	write_taglen($filtered_taglen, $num_filtered, "$outdir/filtered_taglengths.csv");
	#write out soft clipped tag lengths
	write_taglen($softclip_features, $num_softclipped, "$outdir/softclip_taglengths.csv");
	#write out chastity failed tag lengths
	write_taglen($chastity_taglen, $num_chastity, "$outdir/chastity_taglengths.csv");
	#write out individual isoform information
	write_isoform($isoforms, $outdir, $num_mirna);
}
close PROJ;


sub findsamfiles {
	if ($File::Find::name =~ /\.sam$/) {
		push(@samfiles, $File::Find::name);
		print STDERR "\t$File::Find::name\n";
		#check corresponding adapter report is also there at the same time
		my $adapter_file = $File::Find::name;
		$adapter_file =~ s/\.sam$//;
		$adapter_file .= "_adapter.report";
		warn "Adapter report file $adapter_file not found for sam file $File::Find::name, adapter related reports will not be available.\n" unless -e $adapter_file;
	}
}

sub process_file {
	my $file = shift;
	my ($num_reads, $num_unaligned, $num_filtered, $num_overfiltered, $num_softclipped) = (0,0,0,0,0);
	my $num_chastity = 0; #number of aligned, filtered reads that failed chastity
	my %feature_species; #a hash of feature types from annotation priority, each with a hash of species name => read count
	my %crossmapped; #a hash of miRNA species having crossmapped alignments
	my %filtered_taglen; #a hash of tag lengths, each with hashes from annotation priority => read count
	my %softclipped; #a hash of filtered_taglen for softclipped perfect alignments
	my %chastity_taglen; #a hash of filtered_taglen for chastity failed alignments that passed filters
	my %mirna_species; #list of miRNAs seen in sample without regard to XD specifics, used to see miRNA diversity in the sample
	my %isoforms; #list of miRNA isoforms by miRNA name; information contained is coordinate, sequence, read count, crossmapped flag

	open IN, "<$file" or die "Can't open sam file: $file: $!\n";
	while (my $line = <IN>) {
		chomp $line;
		next if $line =~ /^@/;
		
		$num_reads++;
		my ($id, $bitflag, $chr, $start, $mapq, $cigar, $mate, $matepos, $isize, $seq, $qual, $tags) = split(/\t/, $line, 12);
		if ($chr eq '*') {
			$num_unaligned++;
			next;
		}
		
		#tags used: NM, X0, X1, XA, XC, XI, XD
		my ($nm) = $tags =~ /NM:i:(\d+)/;
		my ($x0) = $tags =~ /X0:i:(\d+)/;
		my ($x1) = $tags =~ /X1:i:(\d+)/;
		my ($xa) = $tags =~ /XA:Z:([^\t]*)/;
		my ($xc) = $tags =~ /XC:Z:([^\t]*)/;
		my ($xi) = $tags =~ /XI:Z:([^\t]*)/;
		my ($xd) = $tags =~ /XD:Z:([^\t]*)/;

		die "$file is unannotated. Run the annotation script first" unless $xc ne '';

		next if $nm > $MAX_NM;
		next if $x0 > $MAX_X0;
		$num_filtered++;
		$num_overfiltered++ if $x0 <= $MAX_X0 && $x0 + $x1 > $MAX_XA;
		$num_chastity++ if (($bitflag & 0x0200) == 0x0200);
		
		#split the multi-coordinate tags
		my @xc = split(';', $xc);
		my @xi = split(';', $xi);
		my @xd = split(';', $xd);
		$start = -1 * $start if (($bitflag & 0x0010) == 0x0010); #incorporate strand information into coordinate the same way the XA string does
		my $xa = "$chr,$start,$cigar,$nm;$xa"; #add the coords from the sam line into the xa string to process all coords together
		my @xa = split(';', $xa);

		#first, cycle through checking for miRNA cross-mapping
		my @cross_map_indices = check_crossmap(\@xc, \@xi, \@xd, $x0, \@xa);
		#now grab the "xi,xd" annotation from @xi and @xd using the indices, and sort them
		my @cross_map;
		foreach my $cross_map_index (@cross_map_indices) {
			push(@cross_map, "$xi[$cross_map_index],$xd[$cross_map_index]");
			$mirna_species{$xi[$cross_map_index]} = 1;
		}
		#get the appropriate annotation for the read, which will be used if the read is not crossmapped
		my ($annot_xc, $annot_xi, $annot_xd, $taglength, $softclip, $annot_chr, $annot_pos, $annot_strand) = get_annotation(\@xc, \@xi, \@xd, $x0, \@xa);
		$mirna_species{$annot_xi} = 1 if $annot_xc eq 'miRNA';

		if (scalar @cross_map) {			
			foreach my $annotation (@cross_map) {
				$crossmapped{$annotation} = 0; #for each single miRNA, just add the species name to the hash to get the read counts later
			}
			$crossmapped{join(';', sort @cross_map)}++;
		}

		my $annot_key = scalar @cross_map ? 'crossmapped' : get_annot_key($annot_xc, $annot_xd);
		#dump softclipped alignments to a separate place first
		if ($softclip) {
			$num_softclipped++;
			$softclipped{$taglength}{$annot_key}++;
		}
		else {
			unless (scalar @cross_map || $ANNOTATION_PRIORITY{$annot_key} > 90) { #if @cross_map, then $ANNOTATION_PRIORITY{crossmapped} won't be checked
				$feature_species{$annot_xc}{"$annot_xi,$annot_xd"}++;
			}
			$filtered_taglen{$taglength}{$annot_key}++;
			write_bed($annot_chr, $annot_pos, $taglength, $annot_strand) unless scalar @cross_map; #crossmapped reads will be processed for BED when parsing for isoform information
			$chastity_taglen{$taglength}{$annot_key}++ if (($bitflag & 0x0200) == 0x0200);
		}
		
		#record isoform information
		if ($annot_key eq 'crossmapped') {
			foreach my $i (@cross_map_indices) {
				my ($c_xc, $c_xi, $c_xd, $c_taglength, $c_softclip, $c_chr, $c_pos, $c_strand) =  get_annotation([$xc[$i]], [$xi[$i]], [$xd[$i]], 1, [$xa[$i]]);
				record_isoform(\%isoforms, $c_xi, $c_chr, $c_pos, $c_taglength, $c_strand, $seq, 1, $c_xd);
				write_bed($c_chr, $c_pos, $c_taglength, $c_strand);
			}
		}
		elsif ($ANNOTATION_PRIORITY{$annot_key} < 10) {
			record_isoform(\%isoforms, $annot_xi, $annot_chr, $annot_pos, $taglength, $annot_strand, $seq, 0, $annot_xd);
		}
	}
	close IN;
	#close BED file file handles
	foreach my $chr (keys %BED) {	
		foreach my $strand (keys %{$BED{$chr}}) {
			close $BED{$chr}{$strand};
			system "gzip -f $beddir/$chr\_$strand.txt";
		}
	}
	#clear the BED hash
	%BED = ();

	#if there are any miRNA species in feature_species that are also in crossmapped, take the entry out of feature_species and move it to crossmapped
	foreach my $crossmapped_species (keys %crossmapped) {
		if (exists $feature_species{miRNA}{$crossmapped_species}) {
			$crossmapped{$crossmapped_species} = $feature_species{miRNA}{$crossmapped_species};
			delete $feature_species{miRNA}{$crossmapped_species};			
		}
	}

	return (\%feature_species, \%mirna_species, \%filtered_taglen, \%crossmapped, $num_reads, $num_filtered, $num_unaligned, $num_overfiltered, $num_softclipped, $num_chastity, \%softclipped, \%chastity_taglen, \%isoforms);
}

sub check_crossmap {
	#takes the XC, XI, XD tags for X0 alignments and checks if there are cross-mapped miRNAs
	#if there are none, return empty array
	#if there are cross-maps, return an array of indices corresponding to the input arrays that contain crossmapped entries
	my $xc = shift;
	my $xi = shift;
	my $xd = shift;
	my $x0 = shift;
	my $xa = shift;
	my @crossmaps;
	my @xi_mir; #collect miRNA data with trimmed names (-d nomenclature as described below)
	
	#in case XA is overfiltered, get the smaller of X0 and length of XC array
	$x0 = scalar @{$xc} if scalar @{$xc} < $x0;

	#read only tags corresponding to the X0 alignments
	for (my $i = 0; $i < $x0; $i++) {
		my ($chr, $pos, $cigar, $nm) = split(',', $xa->[$i]);
		next unless $nm <= $MAX_NM;

		if ($xc->[$i] eq 'miRNA') {
			my ($c_xc, $c_xi, $c_xd, $c_taglength, $c_softclip, $c_chr, $c_pos, $c_strand) =  get_annotation([$xc->[$i]], [$xi->[$i]], [$xd->[$i]], 1, [$xa->[$i]]);
			next if $c_softclip;

			my $mirna = join('-', (split('-', $xi->[$i]))[0..2]); #strip off the -d in the a-b-c-d nomenclature in XI, that denotes separate names for the same functional miRNA, since they are not considered cross-mapped
			next if (grep {$_ eq $mirna} @xi_mir);

			push(@xi_mir, $mirna);
			push(@crossmaps, $i);
		}
	}
	#if @crossmaps only has 1 entry, that means that miRNA was not crossmapped;
	return () if scalar @crossmaps <= 1;
	return @crossmaps;
}

sub get_annotation {
	#takes the XC, XI, XD tags for X0 alignments and returns the annotation with the highest priority
	my $xc = shift;
	my $xi = shift;
	my $xd = shift;
	my $x0 = shift;
	my $xa = shift;
	
	#in case XA is overfiltered, get the smaller of X0 and length of XC array
	$x0 = scalar @{$xc} if scalar @{$xc} < $x0;

	my ($xc_ann, $xi_ann, $xd_ann);
	my ($chr, $pos, $cigar, $nm); #fields in the XA string
	my $strand = '+';
	my $taglength; #the length of tag that aligned, after parsing the cigar string
	my $softclip; #0 or 1 indicating whether there was softclipping
	
	my $prev_xtag = "";
	my $index_used = 0;
	#read only tags corresponding to the X0 alignments
	for (my $i = 0; $i < $x0; $i++) {
		my $xtag = get_annot_key($xc->[$i], $xd->[$i]);
		($chr, $pos, $cigar, $nm) = split(',', $xa->[$i]);
		if ($ANNOTATION_PRIORITY{$xtag} < $ANNOTATION_PRIORITY{$prev_xtag} && $nm <= $MAX_NM) {
			$prev_xtag = $xtag;
			$index_used = $i;
		}
	}
	
	$xc_ann = $xc->[$index_used];
	$xi_ann = $xi->[$index_used];
	$xd_ann = $xd->[$index_used];
	($chr, $pos, $cigar, $nm) = split(',', $xa->[$index_used]);
	$pos =~ s/\+//; #strip off the + for positive strand coordinates if necessary
	if ($pos < 0) {
		$strand = '-';
		$pos *= -1;
	}
	$softclip = $cigar =~ /S/;
	$taglength = read_cigar($cigar);

	return ($xc_ann, $xi_ann, $xd_ann, $taglength, $softclip, $chr, $pos, $strand);
}

sub read_cigar {
	#parses cigar string and returns length of sequence
	my $cigar = shift;
	my $taglength = 0;
	
	my @cigar_ops = split(/([MIDNSHP])/, $cigar); #wrapping the regex in () keeps the delimiter as an array element
	while (scalar @cigar_ops) {
		my $numbases = shift @cigar_ops;
		my $op = shift @cigar_ops;
		$taglength += $numbases if $op =~ /[MDNP]/;
	}

	return $taglength;
}

sub get_annot_key {
	#returns the key used in ANNOTATION_PRIORITY
	my $xc = shift;
	my $xd = shift;
	$xd =~ s/,.*//; #remove specific information after the xd tag (mature and star strand miRNAs have the MIMAT ID appended to the end of xd, after a comma)
	return $xd unless $xd eq '';
	return $xc;
}

sub read_adapter {
	my $adapter_file = shift;

	my $adapter_category_index = 0;
	my @adapter_counts;
	#initialize @adapter_counts to 0
	for (my $i = 0; $i < scalar(@ADAPTER_CATEGORIES); $i++) {
		$adapter_counts[$i] = 0;
	}
	#return 0'ed array if there is no adapter_file
	return @adapter_counts unless -e $adapter_file;
	open IN, "<$adapter_file" or die "Can't open adapter report file: $adapter_file: $!\n";
	#adapter reports are ordered by adapter position
	while (<IN>) {
		chomp;
		my ($adapter_pos, $count) = split ' ';
		$adapter_category_index++ if $adapter_pos > $ADAPTER_CATEGORIES[$adapter_category_index] && $ADAPTER_CATEGORIES[$adapter_category_index] >= 0;
		$adapter_counts[$adapter_category_index] += $count;
	}
	close IN;

	return @adapter_counts;
}

sub write_feature_species {
	my $feature_species = shift;
	my $outdir = shift;
	my $num_mirna = shift; #total number of miRNA, for simple coverage normalization

	foreach my $feature (keys %{$feature_species}) {
		my $feature_name = $feature;
		$feature_name =~ s/\s/_/g;

		my @sorted_names = sort {$feature_species->{$feature}{$b} <=> $feature_species->{$feature}{$a}} keys %{$feature_species->{$feature}};
		open FEAT, ">$outdir/$feature_name.txt" or die "Can't write out feature list: $outdir/$feature_name.txt: $!\n";
		foreach my $name (@sorted_names) {
			if ($feature_name eq 'miRNA') {
				my $reads = $feature_species->{$feature}{$name};
				my $pct = sprintf('%.2f%%', $reads * 100 / $num_mirna);
				print FEAT "$name $reads $pct\n";
			}
			else {
				print FEAT "$name $feature_species->{$feature}{$name}\n";
			}
		}
		close FEAT;
	}
}

sub write_mirna_locus {
	my $mirna_locus = shift;
	my $outdir = shift;
	my $num_mirna = shift;

	my @sorted_mirnas = sort {$mirna_locus->{$b} <=> $mirna_locus->{$a}} keys %{$mirna_locus};
	open FH, ">$outdir/mirna_species.txt" or die "Can't write out to file: $outdir/mirna_species.txt: $!\n";
	foreach my $mirna (@sorted_mirnas) {
		my $reads = $mirna_locus->{$mirna};
		my $pct = sprintf('%.2f%%', $reads * 100 / $num_mirna);
		print FH "$mirna $reads $pct\n";
	}
	close FH;
}

sub write_crossmapped {
	my $crossmapped = shift;
	my $outdir = shift;
	my $num_mirna = shift; #total number of miRNA, for simple coverage normalization

	my @sorted_names = sort { $crossmapped->{$b} <=> $crossmapped->{$a} } keys %{$crossmapped};

	open CROSS, ">$outdir/crossmapped.txt" or die "Can't write out crossmapped list: $outdir/crossmapped.txt: $!\n";
	foreach my $name (@sorted_names) {
		next unless $crossmapped->{$name} > 0;
		my $reads = $crossmapped->{$name};
		my $pct = sprintf('%.2f%%', $reads * 100 / $num_mirna);
		print CROSS "$name $reads $pct\n";
	}
}

sub write_taglen {
	my $feature_taglen = shift;
	my $total_filtered = shift;
	my $outfile = shift;

	my @features = ('crossmapped'); #add special annotation type for all crossmapped reads
	foreach my $feature (keys %ANNOTATION_PRIORITY) {
		push(@features, $feature) if $ANNOTATION_PRIORITY{$feature} < 100;
		#order for @features doesn't matter as long as it's consistent throughout the column
		#data will be fed into R which uses column headings to refer to data
	}
	open TAGLEN, ">$outfile" or die "Can't write out taglength stats: $outfile: $!\n";
	print TAGLEN "taglen";
	foreach my $feature (@features) {
		my $name = $feature;
		$name =~ s/\s/_/g; #replace spaces with underscores
		$name =~ s/^(\d)/p$1/; #add a p (prime) in front of 3' and 5' UTR names for easier parsing with R
		print TAGLEN ",$name";
	}
	print TAGLEN "\n";
	my @taglens = sort {$a <=> $b} keys %{$feature_taglen};
	my $mintag = $taglens[0];
	my $maxtag = $taglens[-1];
	for (my $taglen = $mintag; $taglen <= $maxtag; $taglen++) {
		print TAGLEN "$taglen";
		foreach my $feature (@features) {
			my $pct = 0;
			if (exists $feature_taglen->{$taglen} && exists $feature_taglen->{$taglen}{$feature}) {
				$pct = $feature_taglen->{$taglen}{$feature} * 100 / $total_filtered;
			}
			print TAGLEN ",$pct";
		}
		print TAGLEN "\n";
	}
	close TAGLEN;
}

sub write_bed {
	my $chr = shift;
	my $pos = shift;
	my $taglength = shift;
	my $strand = shift;

	unless (exists $BED{$chr}{$strand}) {
		open($BED{$chr}{$strand}, ">$beddir/$chr\_$strand.txt");
	}

	my $bed_start = $pos -1;
	my $bed_end = $bed_start + $taglength; #coordinates comply with BED's 0 indexed, end exclusive numbering
	print {$BED{$chr}{$strand}} "$chr\t$bed_start\t$bed_end\n";
}

sub write_isoform {
	my $isoforms = shift;
	my $outdir = shift;
	my $num_mirna = shift;

	my @sorted_mirnas = sort {$isoforms->{$a} cmp $isoforms->{$b}} keys %{$isoforms};
	open FH, ">$outdir/isoforms.txt" or die "Can't write out to file: $outdir/isoforms.txt: $!\n";
	foreach my $mirna (@sorted_mirnas) {
		foreach my $coord (sort sort_coord_str keys %{$isoforms->{$mirna}}) {
			my ($chr, $start, $end, $strand) = $coord =~ /(.*): (.*) - (.*); (.*)/;
			my $seq = $isoforms->{$mirna}{$coord}{seq};
			my $count = $isoforms->{$mirna}{$coord}{count};
			my $cross = $isoforms->{$mirna}{$coord}{cross};
			my $annot = $isoforms->{$mirna}{$coord}{annot};
			print FH "$mirna $chr $start $end $strand $seq $count $cross $annot\n";
		}
	}
	close FH;
}

sub sort_coord_str {
	my ($a_chr, $a_start, $a_end, $a_strand) = $a =~ /(.*):(.*)-(.*):(.*)/;
	my ($b_chr, $b_start, $b_end, $b_strand) = $b =~ /(.*):(.*)-(.*):(.*)/;
	
	if ($a_chr ne $b_chr) {
		return $a_chr <=> $b_chr;
	}
	elsif ($a_start != $b_start) {
		return $a_start <=> $b_start;
	}
	else {
		return $a_end <=> $b_end;
	}
}

sub record_isoform {
	#adds isoform information about a miRNA read to a hash
	my $isoforms = shift;
	my ($xi, $chr, $pos, $taglength, $strand, $seq, $crossmapped_flag, $xd) = @_;

	my $coord_end = $pos + $taglength - 1;
	my $coord = "$chr: $pos - $coord_end; $strand"; #standard format for coordinate in isoform hash
	if (exists $isoforms->{$xi} && exists $isoforms->{$xi}{$coord}) {
		$isoforms->{$xi}{$coord}{count}++;
		$isoforms->{$xi}{$coord}{cross} = 1 if $crossmapped_flag == 1;
		$isoforms->{$xi}{$coord}{annot} = $xd;
	}
	else {
		$isoforms->{$xi}{$coord}{seq} = $seq; 
		$isoforms->{$xi}{$coord}{count} = 1; 
		$isoforms->{$xi}{$coord}{cross} = $crossmapped_flag; 
		$isoforms->{$xi}{$coord}{annot} = $xd; 
	}
}
