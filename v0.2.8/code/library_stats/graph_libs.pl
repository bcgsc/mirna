#!/usr/bin/env perl
use strict;
use Getopt::Std;
use Pod::Usage;
use File::Find;
use File::Basename;
use File::Copy;

use vars qw($opt_p);
getopts("p:");

my $usage = "$0 -p project_directory\n";
die "$usage" unless $opt_p;

my $r = get_config();

print STDERR "Searching for .sam files to in project directory...\n";
my @samfiles;
find(\&findsamfiles, $opt_p);
print STDERR "Done\n";
die "No sam files found in $opt_p" unless scalar @samfiles;

my $dir = dirname(__FILE__); #R scripts are in the same directory as this script

mkdir "$opt_p/graphs" unless -e "$opt_p/graphs";
mkdir "$opt_p/graphs/tags" unless -e "$opt_p/graphs/tags";
mkdir "$opt_p/graphs/adapter" unless -e "$opt_p/graphs/adapter";

foreach my $samfile (@samfiles) {
	my $samdir = dirname($samfile);
	my $filename = basename($samfile);
	$filename =~ s/\.[bs]am$//;
	my ($lib, $index) = split('_', $filename);
	$index = '' unless defined $index;

	my $datadir = "$samdir/$filename\_features";
	my $filtered_file = "$datadir/filtered_taglengths.csv";
	my $softclip_file = "$datadir/softclip_taglengths.csv";
	my $chastity_file = "$datadir/chastity_taglengths.csv";
	my $adapter_file = "$samdir/$filename\_adapter.report";
	
	system "$r $dir/taglengths.R $datadir $filename $filtered_file tags";
	system "$r $dir/taglengths.R $datadir $filename $softclip_file softclip";
	system "$r $dir/taglengths.R $datadir $filename $chastity_file chastity";

	system "$r $dir/adapter.R $datadir $filename $adapter_file adapter" if -e $adapter_file;

	#make a copy of tags and adapter under the graphs directory for each access
	copy("$datadir/$filename\_tags.jpg", "$opt_p/graphs/tags/$filename\_tags.jpg");
	copy("$datadir/$filename\_adapter.jpg", "$opt_p/graphs/adapter/$filename\_adapter.jpg");
}

my $proj = basename($opt_p); #use project name for saturation graph title
my $saturation_source = "$opt_p/alignment_stats.csv";
system "$r $dir/saturation.R $opt_p/graphs/$proj\_saturation.jpg ".uc($proj)." $saturation_source";

sub findsamfiles {
	if ($File::Find::name =~ /\.sam$/ || $File::Find::name =~ /\.bam$/) {
		#skip specialized analyses within _features directories
		next if $File::Find::name =~ /_features/;
		#skip files in obsoleted directory
		next if $File::Find::name =~ /obsoleted/;
		push(@samfiles, $File::Find::name);
		print STDERR "\t$File::Find::name\n";
	}
}

sub get_config {
	my $dir = dirname(__FILE__);
	my $config_file = "$dir/../../config/pipeline_params.cfg";
	
	open CONFIG, $config_file or die "Could not find config file in default location ($config_file)";
	my @config = <CONFIG>;
	close CONFIG;
	chomp @config;
	
	my ($r) = [grep(/^\s*Rscript/, @config)]->[0] =~ /^\s*Rscript\s*=\s*(.+)/ or die "No entry found for Rscript binary in config file.";
	return $r;
}
