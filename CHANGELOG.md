# 1.0.0
## Changed
 - Refactored and Rewritten the entire pipeline in (pure) Python
 - Added a column to isoforms.txt for isomir nomenclature
## Fixed
 - Fixed a bug in the annotation module where it assumed that the coding region spans all exons, and so over-reported the number of 5' and 3' UTR exons, and under-reported other features.
 - Fixed a query bug that caused the pipeline to miss some of the snoRNA genes in the UCSC database, and so under-report snoRNA hits.
 - Fixed a typo that caused the pipeline to miss all Low-complexity repeat hits
## Added
 - Added functionality to retrieve miRNA GFF3 files directly from miRBase FTP site, and thereby removed dependency on local mirbase database
 - Added Dockerfile for containerization of the pipeline. The container is pushed to DockerHub.
 - Added Nextflow workflow which runs the pipeline container pulled from DockerHub with Singularity

# v.0.2.8 
Changed read alignment end coordinates in the isoform.quantification.txt file to be inclusive. They had been exclusive (i.e. offset by 1).

# v.0.2.7
Add compatibility with miRBase v20 schema updates.
Update expression matrix generation scripts to allow samples with 0 miRNA reads by writing a warning to STDERR instead of stopping script.
Add script to write expression matrix of miRNA mature strands.
Add script to write expression matrix of miRNA mature strands using TCGA data (Level 3 archive and ADF reference).
Moved adapter trimming preprocessing to a separate location (http://www.bcgsc.ca/platform/bioinfo/software/adapter-trimming-for-small-rna-sequencing), removed extraneous scripts.
Update graphing calls to use Rscript instead of R in batch mode.
Update saturation graph to draw line of best fit only when there are enough samples to do so.

# v.0.2.6
Create new graphing script to show miRNA read count vs species diversity, using all samples sequenced, including failed and merged repool samples.
Place a copy of the graphs under the main project directory as well as inside the sample _features directory for easier access.

Cross-mapped reads are now used for expression matrix and enumerating miRNA species diversity; they are double-counted.
This means that tcga.pl, which creates the DCC quantification files, must use miRNA.txt instead of mirna_species.txt.

TCGA formatted data files have an updated format:
	headers are in the files
	all miRNAs are listed in the miRNA quanitification file, rather than just those miRNAs with non-zero expression, so they are all the same length now (isoform quantification file lengths are unchanged - they are dependent on the reads in the sample)
	normalized read counts use floating point numbers (6 decimal places) rather than rounding off to integer
TCGA ADF reports miRNA accession from miRBase as well.

Add a new custom output, for BCGSC, which helps track how samples are processed within the BCGSC.

Bug fix:
TCGA ADF file now reports MIMAT accessions separately instead of grouping it with coordinates. This means that the alternate mature strand will always take up 2 columns and will not cause column misalignment.

Configuration changes:
Path to samtools binary used is now referenced by config/pipeline_params.cfg.
Use 64bit R instead of 32bit in config/pipeline_params.cfg.
GRCh37 reference databases added to config/db_connections.cfg.

Minor script improvements:

expression_matrix.pl:
Columns in expression matrix are now sorted by sample ID
Add "Gene" as the header for the first column instead of leaving it blank.
Only include miRNAs in the expression matrix where there is coordinate and mature strand information in miRBase, instead of having an extra row where all values will be 0.
Normalized and log counts are rounded to 6 decimal places.

alignment_stats.pl:
Update sorting order of isoforms result file such that sorting is done alphanumerically instead of numerically.
Fill in 0 instead of blank for empty values.
Instead of failing when an adapter report is not found, use 0 for adapter counts and print a warning.
When searching for input files, specifically ignore anything where the pathname includes "obsoleted".

graph_libs.pl
When searching for input files, ignore anything where pathname includes "obsoleted", and ignore files under the _features directories, which are specialized analyses that doesn't need graphs.
Search for .bams as well as .sams for indications of a sample.


# v0.2.5
Bug fix in incorporating MIMAT accessions into star strands of miRNAs. MIMAT got incorporated into the XI tag leading to no counts of any star strand in alignment_stats.csv.

Bug fix in reannotating a sam file where previous annotation was not stripped off if the annotation had white space. Now the replacement conforms to the correct syntax where only tabs are not allowed in the annotation.

# v0.2.4
Major bug fix in annotation of mature and star strands of miRNAs on the - strand. On the - strand, the genomic coordinates of the mature and star strands are equal to the precursor end coordinate minus the miRBase mature_from and mature_to offsets, rather than adding them to the precursor start coordinate as with + strand miRNAs.

Bug fix in counting cross-mapped miRNAs.

Incorporate MIMAT accessions with mature and star annotations.

Generate 1 set of BED files for each sample instead of for an entire project, so coverage files can be generated for individual samples.

Add extra column to alignment_stats.csv - adapter dimer %.

# v0.2.3
Added a new stats script to generate a miRNA expression matrix from the sample data.

Folded code from annotate.pl to module Annotate.pm so that all "generic annotation code" is in the module. This improves code reuse in a new script to annotate .peaks files generated by FindPeaks (new feature in future release).

Use the Perl File::Copy module's move function to overwrite the unannotated files with the annotated versions, instead of making a system call to the shell mv. This is an attempt to fix a not-as-yet consistently reproducible bug where the annotated file remains with its temporary name.

Updated graph such that the increased reporting categories accumulated from previous versions does not overwhelm the colour palette. Detailed category data is now only available in alignment_stats.csv. The graphs summarize the groups down such that they can be visually separated by eye.

alignment_stats.pl:
- bug fix for counting cross-mapped miRNAs. Fix makes the script only look at BWA's optimal matches (X0) and not X1 when considering cross-mapping.
- write out isoform information about where along the miRNA different reads are found.
- minor bug fix; check each annotation's coordinate string to make sure that each NM field does not exceed max number of allowed mismatches MAX_NM

Custom output:
A new set of scripts takes the output from alignment_stats.pl and formats it to conform to standards dictated by various projects. Currently only applies to TCGA (http://tcga.cancer.gov/) miRNA data.

# v0.2.2
Initial public release. Added a README.TXT at the base level directory that explains how to run the pipeline.

Add new reporting features to analysis.
- instead of just annotating RNAs in RepeatMasker, annotate every repeat class, and split each into a separate category
- report the miRNA species represented without regard to specific annotation
- enumerate the number of miRNA species in the sample with a >= 10x depth
See code/library_stats/README.TXT for specific reporting details.

Instead of splitting up BED files generated by category, leave them all together so all peaks will show up in 1 wig. In a future release, the peaks themselves will be re-annotated instead of relying on the annotations of the reads which make them up.

Minor improvement in coding style for referencing modules.
 - instead of looking for the module in ".", which requires the code to be run from the script directory, the script only requires the module to be in the same directory the script is in, so now the script can be called from any path.

# v0.2.1
Begin integration of novel miRNA prediction.
alignment_stats.pl now generates a set of bed files from the read coordinates for use in peak detection by FindPeaks. (http://sourceforge.net/apps/mediawiki/vancouvershortr/index.php?title=Main_Page#FindPeaks_4.0)
These peaks will be used for predicting novel miRNAs.

Minor bug fix for reading the additional coordinate information in the XA tag.

# v0.2.0
Increase reporting categories.
- annotated tRNAs separately from other RNAs in RepeatMasker
- add meta-stat to check the expression profile of chastity failed reads
- report number of reads with 3' adapter in positions 0, 1-14, 15-25, 26-35, and 35+

Update meta-stat regarding overfiltered (X0+X1 > samse -n) stat to reflect up to 10 alignments reported by LIMS
Fixed minor bug where the same stat was reported twice in the summary csv.

Modify colour scheme of graphs to improve readability with increased categories
Add readme to document the annotation output.
Add extra graphs to report 3' adapter location along read.
Add a normalization to miRNA expression files (miRNA.txt and crossmapped.txt) that shows read count as a percentage of all miRNA found in the library.

# v0.1.0
First production version of pipeline for annotating miRNA libraries.
