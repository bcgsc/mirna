DESCRIPTION OF DATA FILES FROM MIRNA PROFILING

I. alignment_stats.csv
II. graphs
III. Sample specific data files
IV. TCGA formatted output
V. Expression matrix

I. alignment_stats.csv
alignment_stats.csv in the base directory is a summary of results for each sample. The individual fields are described as follows.

Library, Index: Identifier of the sample.
Total Reads: number of all reads used in alignment that can be associated with the sample. In the case of indexed libraries, this means the set of reads where the first N bases map back to the index sequence of length N. Note that since reads with less than 15bp of biological sequence were not sent for alignment, this number does not reflect those reads.
% Adapter dimers: number of adapter dimers as percentage of all reads associated with library, including reads with adapter before the 15th base which were not includedin the "Total Reads" column.
Adapter dimers: number of reads where the 3' adapter was found at the beginning of the read, implying 5' and 3' adapters formed a dimer with no biological RNA in the sequence. Discarded without alignment.
Adapter at 1-14bp: number of reads where the 3' adapter was found at the first 14bp of the read, meaning that there is less 15bp of biological sequence. Discarded without alignment.
Adapter at 15-25bp; Adapter at 26-35bp, Adapter after 35bp: number of reads with 3' adapter found after n bp of sequence.
Aligned Reads Post-Filter: number of reads that are considered for analysis. This is the number of reads that align perfectly to the genome, and map to 3 or fewer places. This number will be the basis of percentages for the feature counts.
Unaligned Reads: number of reads that could not be aligned to the genome at all. Total Reads - (Aligned Reads Post-Filter + Unaligned Reads) = number of reads that aligned but with any combination of mismatches, indels, mapped to too many regions, etc.
% Aligned Reads, %Unaligned Reads: Aligned and Unaligned Reads expressed as a percentage of Total Reads.
Filtered Reads without XA: ***Used for internal development*** Number of reads where extra mapped positions could not be used because they were not reported by bwa when X0 + X1 > N in `bwa samse -n N`. This is a subset of the Aligned Reads Post-Filter.
Softclipped Reads: ***Used for internal development*** Number of reads that were perfectly aligned to the genome after bwa softclips part of the sequence. These alignments ARE NOT used in the core reports of reads mapped to features.
Chastity Failed Reads Post-Filter: ***Used for internal development*** Number of reads that aligned perfectly and passed all filters, but failed chastity check. These alignments ARE used in the core reports of reads mapped to features.
miRNA Species: number of miRNA species covered by at least 1 read.
miRNA Species Covered by >= 10 Reads: number of miRNA species covered by at least 10 reads. Crossmapped reads will be counted once for each miRNA it aligns to (ie. 1 read will be counted as coverage for multiple miRNAs).
Total miRNA: number of reads aligning to miRNAs, this number is the sum of Crossmapped miRNA, mature miRNA, * miRNA, precursor miRNA, miRNA loop, and unannotated miRNA
Crossmapped miRNA: number of reads that align to crossmapped miRNAs, regardless of whether the read itself is crossmapped. A crossmapped miRNA has reads that perfectly to more than 1 miRNA species, regardless of location within the miRNA. These reads are not counted in the mature, *, precursor, loop and unannotated miRNA numbers. A separate sum of reads that are and are not crossmapped can be found by using the sample specific data file crossmapped.txt, described below.
mature miRNA: number of reads mapped to the mature strand of miRNAs
* miRNA: number of reads mapped to the * strand of miRNAs
precursor miRNA: number of reads that either 1. do not map to the mature strand, but map to the region between the end of the mature strand and the end of the precursor miRNA, or 2. do not map to the * star strand, but map to the region between the end of the * strand and the end of the precursor miRNA
miRNA loop: number of reads that do not align to the mature strand, but align to the 5 bases after the mature strand, where the loop between the mature strand and * strand would be.
unannotated miRNA: number of reads that align to the region of the precursor miRNA that would appear to be on the opposite side of the mature strand, but miRBase does not have the * strand annotation. Therefore these reads cannot be annotated as * strand since there is no information regarding which part of the precursor miRNA would be the * strand.
snoRNA: number of reads that map to known snoRNAs. The list of snoRNAs is extracted from UCSC knownGenes and kgXref where geneSymbol is in the form "SNOR*" or the description contains the string " small nucleolar RNA".
tRNA; rRNA; snRNA; scRNA; srpRNA: number of reads that map to regions as listed by UCSC Genome Browser's RepeatMasker with the given repeat_class
Other RepeatMasker RNAs: number of reads that map to regions as listed by UCSC Genome Browser's RepeatMasker with the repeat_class "RNA"
RNA (No CDS): number of reads that map to genes in UCSC Genome Browser's knownGenes that have 0 length CDS regions (excluding snoRNAs matched by the filter described above)
3' UTR, 5' UTR: number of reads that map to 3' and 5' UTR Exons according to UCSC Genome Browser's knownGenes tables, after strand, coordinate, and CDS region has been calculated
Coding Exon, Intron: number of reads that map to the exon and intron regions of UCSC knownGenes
LINE, SINE, LTR, Satellite, RepeatMasker DNA, RepeatMasker Low complexity, RepeatMasker Simple repeat, RepeatMasker Other, RepeatMasker Unknown: number of reads mapping to other repeat_class categories of UCSC RepeatMasker
Unknown: number of reads that map to no features in the above reference datasets
% Total miRNA to Unknown: number of reads corresponding to above features, expressed as a percentage of Aligned Reads Post-Filter

II. graphs
A _saturation.jpg graph in the graphs subdirectory shows the number of reads aligned to miRNAs, and the number of miRNAs covered by these reads. In the subdirectories there are 2 types of graphs available. tags graphs the percentage of reads aligning to various types of non-coding RNAs, and the read lengths of these aligned reads. adapter graphs the read length distribution of all reads, including adapter-adapter dimers with no RNA.

III. Sample specific data files
For each of the features summarized in alignment_stats.csv, files were created which lists the specific gene name expressed, and the expression level as a read count, sorted by expression.
Coding_exon.txt contains the list of all exons originating from UCSC knownGene. The "RNA (No CDS)" set is differentiated from other exons with the "No CDS" annotation.
crossmapped.txt contains all the miRNAs that have been crossmapped by any reads. This set includes the count of reads that did crossmap, and the set of reads that did not crossmap themselves, but map to reads that have been crossmapped. These 2 categories are distinguished by the line either containing 1 miRNA, or mulitple comma separated miRNA names.
miRNA.txt contains the list of all uniquely mapped miRNAs, with an extra annotation denoting the region of the miRNA the read mapped to.
miRNA expression files, miRNA.txt and crossmapped.txt have a 3rd column which is the read count as a percentage of "Total miRNA" from alignment_stats.csv.
mirna_species.txt is a summary of the miRNA.txt file, with specific annotations about miRNA region left out. This means that each count is the number of reads aligning anywhere along the miRNA. Crossmapped reads will be counted once for each miRNA it aligns to (ie. 1 read will be counted as coverage for multiple miRNAs).
isoforms.txt includes miRNAs in the miRNA.txt and crossmapped.txt files, split by individual reads that make up the miRNA. Column definitions are: miRNA, chr, start, end, strand, read sequence, read count, crossmapped (1 = yes, 0 = no), annotation.
A bed directory is created for each sample. These directories contain text files in BED format showing coverage of all reads which pass the alignment filters. The text files can be used to generate wig coverage files. They are also sorted and filtered so that they can be used for peak discovery and novel miRNA prediction.
All other .csv files are used for generating graphs.

IV. TCGA formatted output
Expression of miRNA and isoforms are also reported in a format compatible with that used by The Cancer Genome Atlas. These results are in the _features/tcga subdirectory. The specifications of this format are available at https://wiki.nci.nih.gov/display/TCGA/miRNASeq

V. Expression matrix
A tab separated expression matrix of samples by miRNAs are available in the base directory. These are expn_matrix.txt, expn_matrix_norm.txt, expn_matrix_norm_log.txt.
expn_matrix.txt is a raw count of the number of reads aligned to a miRNA for each sample, 0 represents 0 reads aligned to a particular miRNA.
expn_matrix_norm.txt normalizes the raw read counts for the different number of reads sequenced for each library such that the sum of normalized expression for any library is always 1 million.
expn_matrix_norm_log.txt a log2 of each value of the normalized file. 0s in the normalized file remain as 0 in the log file as a convenience to allow for numerical processing by programs that don't recognize the string "-Inf". The 0s are actually negative infinity.
Cross-mapped reads are counted once for each miRNA they align to.
