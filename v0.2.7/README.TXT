GSC MiRNA Analysis Pipeline v0.2.7

Aim:
Given a set of aligned reads in 1 or more .sam files, produce an annotated version of the .sam where each read is given an annotation based on its coordinate. Additional summary information about the content of each sample is also generated, including miRNA species and other genomic features found.

To Use:
The code package consists of 3 subdirectories,
apps
code
config.

apps contain external applications required to run the pipeline. Although applications do not have to be stored in this directory, the following binaries must be available on your system: Perl, R. In addition, perl requires the MySQL DBI library. R is used to generate summary graphs, and may be disregarded if graphs aren't desired.

config contains configuration files which points the code to all necessary inputs.
db_connections.cfg contains parameters to access MySQL databases containing UCSC and miRBase information. The db_name is used when providing the database source for various script parameters. You must have a database connection to a miRBase instance and a UCSC database instance for annotations of miRNAs and other non-coding RNAs respectively. The server name of the database is the <host> parameter, and the login and password are the <user> and <password> parameters.
profile.sh points to the appropriate perl to use for the scripts. Run source on this to generate the environment.

code contains all the scripts to run the pipeline. This consists of 2 parts, annotation, and gather library statistics.

First, set up the directory structure such that the scripts can automatically find the files to process. Under the project base directory ({$PROJDIR}), place all your .sam files. These should be named LIBRARY.sam or LIBRARY_INDEX.sam for multiplexed runs. Accompanying each file should be a LIBRARY_adapter.report (or LIBRARY_INDEX_adapter.report) file which is simply a space delimited file which summarizes the adapter trimming done on the reads before alignment. The format is
"read length (from 0 to length of full read)" "number of reads"
eg.
0 1200000
1 123
...
35 90000
36 900000

After the directory structure is set up, run code/annotate/annotate.pl as described by code/annotate/annotate.pl/HOWTO.txt. While the script is running, it will create a version of the .sam named .sam.annot which adds the tags XC, XI, and XD tags to each line, leaving everything else identical. Upon completion of each file, the .sam.annot will overwrite the .sam. NOTE: the bitmap flag in the .sam needs to be in hexidecimal format for bit calculations, the read will not be processed correctly if string format is used.

After annotation, summary stats can be generated. Run the alignment_stats.pl script
> alignment_stats.pl -p {$PROJDIR}
The output of this script is explained in code/library_stats/README.TXT.
Graphs can be generated to provide a visual summary of sample content. Make sure that if running this script through an ssh tunnel, X11 forwarding is enabled so that the drawing libraries in R can be used. eg. > ssh -X $myserver
Run the script
> graph_libs.pl -p {$PROJDIR}
An expression matrix can be built that summarizes the expression of each miRNA gene in each sample in the project. Run the script
> expression_matrix.pl -m miRBase_db_name -o miRBase_species_name -p project_dir
in the same way as alignment_stats.pl.
The equivalent expression matrix that reports read counts for miRNA mature strands rather than miRNA precursors can be generated using
> expression_matrix_mimat.pl -m miRBase_db_name -o miRBase_species_name -p project_dir

For further questions and support requests, contact Andy Chu at achu@bcgsc.ca

Appendix.

Chromosome/Coordinate Labelling
When annotating a file, the alignment coordinates are checked against reference databases for overlap. To account for differences in chromosome naming convention (eg. chr1 vs 1), a "sample" is taken from one of the files to annotate and labelling format is checked.
