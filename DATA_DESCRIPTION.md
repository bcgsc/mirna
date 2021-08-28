# miRNA Profiling Container Output

## I. miRNA Expression

### miRNA.txt

 lists all uniquely mapped miRNAs, with an extra annotation denoting
 the region of the miRNA the read mapped to.

<table>
<thead>
<td><strong>Column</strong></td>
<td>
<strong>Ttile</strong>
</td>
<td><strong>Type</strong></td>
<td>
<strong>Example</strong>
</td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>miRNA feature</td>
<td>
str
</td>
<td>hsa-mir-103a-2,mature,MIMAT0000101</td>
<td>
miRNA feature name
</td>
</tr>
<tr>
<td>
2
</td>
<td>
read count
</td>
<td></td>
<td>
62
</td>
<td>
Number of aligned reads mapped to the miRNA feature
</td>
</tr>
<tr>
<td>
3
</td>
<td>% read count</td>
<td></td>
<td>
0.01%
</td>
<td>Percentage of reads mapped to miRNA that are mapped to the miRNA feature</td>
</tr>
</tbody>
</table>

### crossmapped.txt

 contains all miRNAs that have been crossmapped by any reads. This set
 includes the count of reads that do crossmap, and the set of reads
 that do not crossmap themselves, but map to miRNAs that have been
 crossmapped. These 2 categories are distinguished by the line either
 containing 1 miRNA, or multiple comma-separated miRNA names.

<table>
<thead>
<td><strong>Column</strong></td>
<td>
<strong>Title</strong>
</td>
<td><strong>Type</strong></td>
<td>
<strong>Example</strong>
</td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>
miRNA
feature(s)
</td>
<td>
str
</td>
<td>
hsa-mir-103a-1,mature,MIMAT0000101
or
hsa-mir-103a-1,mature,MIMAT0000101;hsa-mir-107, mature,MIMAT0000104
</td>
<td>
A semicolon-separated list of miRNA feature id(s) that are crossmapped to reads
</td>
</tr>
<tr>
<td>
2
</td>
<td>
read count
</td>
<td>
int
</td>
<td>
55586
</td>
<td>
Number of aligned reads mapped to the crossmapped miRNA feature (s)
</td>
</tr>
<tr>
<td>
3
</td>
<td>
% read
count
</td>
<td></td>
<td>
4.95%
</td>
<td>
Percentage of reads mapped to miRNA that are mapped to the crossmapped miRNA feature(s)
</td>
</tr>
</tbody>
</table>

### mirna_species.txt

 List all mapped miRNA primary transcripts. Each read count is the
 number of reads aligning anywhere along the miRNA. Crossmapped reads
 will be counted once for each miRNA it aligns to

<table>
<thead>
<td><strong>Column</strong></td>
<td>
<strong>Title</strong>
</td>
<td><strong>Type</strong></td>
<td>
<strong>Example</strong>
</td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>miRNA species</td>
<td>
str
</td>
<td>hsa-mir-103a-1</td>
<td>
miRNA species, i.e., miRNA primary transcript, name
</td>
</tr>
<tr>
<td>
2
</td>
<td>
read count
</td>
<td>
int
</td>
<td>
55694
</td>
<td>
Number of aligned reads mapped to the miRNA species
</td>
</tr>
<tr>
<td>
3
</td>
<td>
% read count
</td>
<td></td>
<td>
4.96%
</td>
<td>Percentage of reads mapped to miRNA that are mapped to the miRNA species</td>
</tr>
</tbody>
</table>

### isoforms.txt

 contains miRNAs in miRNA.txt and crossmapped.txt files, split by
 individual reads that make up the miRNA.

<table>
<thead>
<td><strong>Column</strong></td>
<td>
<strong>Title</strong>
</td>
<td><strong>Type</strong></td>
<td>
<strong>Example</strong>
</td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>miRNA_ID</td>
<td>
str
</td>
<td>
hsa-let-7a-1
</td>
<td>
miRNA primary transcript id
</td>
</tr>
<tr>
<td>
2
</td>
<td>
chr
</td>
<td>
str
</td>
<td>
9
</td>
<td>
chromosome name
</td>
</tr>
<tr>
<td>
3
</td>
<td>
start
</td>
<td>
int
</td>
<td>
96938244
</td>
<td>
start chromosomal coordinate
</td>
</tr>
<tr>
<td>
4
</td>
<td>
end
</td>
<td>
int
</td>
<td>
96938267
</td>
<td>
end chromosomal coordinate
</td>
</tr>
<tr>
<td>
5
</td>
<td>
strand
</td>
<td>
str
</td>
<td>
+
</td>
<td>
strand
</td>
</tr>
<tr>
<td>
6
</td>
<td>
isomir
sequence</td>
<td>
str
</td>
<td>
TGAGGTAGTAG GTTGTATAGTTTT
</td>
<td>
isomir sequence
</td>
</tr>
<tr>
<td>
7
</td>
<td>
read
count
</td>
<td></td>
<td>
21
</td>
<td>
read count
</td>
</tr>
<tr>
<td>
8
</td>
<td>
crossmapped
</td>
<td>
int
</td>
<td>
0
</td>
<td>
1 if feature is crossmapped, 0 otherwise
</td>
</tr>
<tr>
<td>
9
</td>
<td>
miRNA_feature
</td>
<td>
str
</td>
<td>
mature,
MIMAT0000062
</td>
<td>
Type of miRNA feature
</td>
</tr>
<tr>
<td>
10
</td>
<td>
isomir
nomenclat ure
</td>
<td>
str
</td>
<td>hsa-let-7a-5p|0|+2|</td>
<td>
mature miRNA id followed by the start and end positions of the isomir relative to those of the reference mature strand, - means isomir starts or ends upstream of the reference mature strand, + downstream.
</td>
</tr></tbody>
</table>

### TCGA-formatted miRNA Expression

 Expression of miRNA and isoforms are also reported in a format
 compatible with that used by The Cancer Genome Atlas (TCGA). These
 results are in the \_features/tcga subdirectory. The specifications of
 this format are available at
 https://wiki.nci.nih.gov/display/TCGA/miRNASeq

 tcga/mirnas.txt

<table>
<thead>
<td><strong>Column</strong></td>
<td>
<strong>Title</strong>
</td>
<td><strong>Type</strong></td>
<td><strong>Example</strong></td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>
miRNA_ID
</td>
<td>
str
</td>
<td>
hsa-let-7d
</td>
<td>
valid miRBASE id
</td>
</tr>
<tr>
<td>
2
</td>
<td>
read_count
</td>
<td>
int
</td>
<td>
1285
</td>
<td>
raw read count
</td>
</tr>
<tr>
<td>
3
</td>
<td>reads_per_million_miRNA_mapped</td>
<td>
float
</td>
<td>314.747706</td>
<td>
normalized expression value
</td>
</tr>
<tr>
<td>
4
</td>
<td>
cross_mapped
</td>
<td>
str
</td>
<td>
N
</td>
<td>Is this feature cross-classified as other miRNA forms</td>
</tr>
</tbody>
</table>

 tcga/isoforms.txt

<table>
<thead>
<td><strong>Column</strong></td>
<td>
<strong>Title</strong>
</td>
<td><strong>Type</strong></td>
<td>
<strong>Example</strong>
</td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>
miRNA_ID
</td>
<td>
str
</td>
<td>
hsa-let-7a-1
</td>
<td>
valid miRBASE id
</td>
</tr>
<tr>
<td>
2
</td>
<td>
isoform_coords
</td>
<td>
str
</td>
<td>hg19:9:96938244-96938267:+</td>
<td>
genomic coordinates of the isoform
</td>
</tr>
<tr>
<td>
3
</td>
<td>
read_count
</td>
<td>
int
</td>
<td>
474
</td>
<td>
read_count
</td>
</tr>
<tr>
<td>
4
</td>
<td>reads_per_million_miRNA_mapped</td>
<td>
float
</td>
<td>
116.101488
</td>
<td>
normalized expression value
</td>
</tr>
<tr>
<td>
5
</td>
<td>
cross_mapped
</td>
<td>
str
</td>
<td>
N
</td>
<td>Is this feature cross-classified as other miRNA forms</td>
</tr>
<tr>
<td>
6
</td>
<td>
miRNA_feature
</td>
<td>
str
</td>
<td>
mature,MIMAT0000062
</td>
<td>
Type of miRNA feature
</td>
</tr>
</tbody>
</table>

## II. Alignment Stat

### alignment_stats.csv

 Summary stat of reads and feature expression in the sample(s)

<table>
<thead>
<td><strong>Column</strong></td>
<td>
<strong>Title</strong>
</td>
<td><strong>Type</strong></td>
<td><strong>Example</strong></td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>
Library
</td>
<td>
str
</td>
<td>
m19542
</td>
<td>
Pool or Library name
</td>
</tr>
<tr>
<td>
2
</td>
<td>
Total Reads &gt;= 15bp
</td>
<td>
int
</td>
<td>
13058648
</td>
<td>
Total number of reads (&gt;= 15bp long). 15 was chosen for BCGSC bams because it is the shortest mature miRNA in miRBase 16, but it may be different for bams from other sources.
</td>
</tr>
<tr>
<td>
3
</td>
<td>
% Adapter dimers
</td>
<td></td>
<td>
9.27%
</td>
<td>
Number of adapter dimers as percentage of all reads associated with library, including reads &lt;15bp long
</td>
</tr>
<tr>
<td>
4
</td>
<td>
Adapter dimers
</td>
<td>
int
</td>
<td>
1791787
</td>
<td>
Number of unaligned reads where 5' and 3' adapters formed a dimer with no biological RNA sequence.
</td>
</tr>
<tr>
<td>
5
</td>
<td>
Adapter at 1-14bp
</td>
<td>
int
</td>
<td>
3911832
</td>
<td>
Number of unaligned reads where there is less than 15bp of biological sequence.
</td>
</tr>
<tr>
<td>
6
</td>
<td>
Adapter at 15-25bp
</td>
<td>
int
</td>
<td>
7156962
</td>
<td>
Number of reads with 3' adapter found at 15-25bp of sequence.
</td>
</tr>
<tr>
<td>
7
</td>
<td>
Adapter at 26-35bp
</td>
<td>
int
</td>
<td>
5636422
</td>
<td>
Number of reads with 3' adapter found at 26-35bp of sequence.
</td>
</tr>
<tr>
<td>
8
</td>
<td>
Adapter after 35bp
</td>
<td>
int
</td>
<td>
825737
</td>
<td>
Number of reads with 3' adapter found after 35bp of sequence.
</td>
</tr>
<tr>
<td>
9
</td>
<td>
Aligned Reads Post Filter
</td>
<td>
int
</td>
<td>
5178136
</td>
<td>
Number of aligned and filtered reads, i.e. they align perfectly to the genome, AND map to 3 or fewer places. This number is the basis of % feature read count calculations.
</td>
</tr>
<tr>
<td>
10
</td>
<td>
% Aligned Reads
</td>
<td></td>
<td>
39.65%
</td>
<td>
Percentage of total reads aligned and filtered.
</td>
</tr>
<tr>
<td>
11
</td>
<td>
Unaligned Reads
</td>
<td>
int
</td>
<td>
600996
</td>
<td>
Number of unaligned reads
</td>
</tr>
<tr>
<td>
12
</td>
<td>
% Unaligned Reads
</td>
<td></td>
<td>
4.60%
</td>
<td>
Percentage of total reads unaligned
</td>
</tr>
<tr>
<td>
13
</td>
<td>
Filtered reads without XA
</td>
<td>
int
</td>
<td>
283218
</td>
<td>
Number of aligned reads with alignments not reported in XA, i.e., they were not reported by bwa when X0 + X1 &gt; N in `bwa samse -n N`.
</td>
</tr>
<tr>
<td>
14
</td>
<td>
Softclipped Reads
</td>
<td>
int
</td>
<td>
0
</td>
<td>
Number of aligned reads that are softclipped, i.e. they align perfectly to the genome after bwa softclips part of the sequence. These alignments ARE NOT counted in the reports of reads mapped to features.
</td>
</tr>
<tr>
<td>
15
</td>
<td>
Chastity Failed Reads Post-Filter
</td>
<td>
int
</td>
<td>
0
</td>
<td>
Number of aligned reads that failed chastity, i.e. they align perfectly and pass all filters, but fail chastity check. These alignments ARE counted in the reports of reads mapped to features.
</td>
</tr>
<tr>
<td>
16
</td>
<td>
miRNA Species
</td>
<td>
int
</td>
<td>
859
</td>
<td>
Number of miRNA species mapped to aligned reads
</td>
</tr>
<tr>
<td>
17
</td>
<td>
miRNA Species
Covered by &gt;= 10
Reads
</td>
<td>
int
</td>
<td>
378
</td>
<td>
Number of miRNA species mapped to &gt;= 10 aligned reads
</td>
</tr>
<tr>
<td>
18
</td>
<td>
Total miRNA
</td>
<td>
int
</td>
<td>
2107670
</td>
<td>
Number of aligned reads mapped to miRNA species. This number is the basis of % read count calculations in miRNA expression files.
</td>
</tr>
<tr>
<td>
19
</td>
<td>Crossmapped miRNA</td>
<td>
int
</td>
<td>
300394
</td>
<td>
Number of aligned reads crossmapped to miRNA species
</td>
</tr>
<tr>
<td>
20
</td>
<td>
mature miRNA
</td>
<td>
int
</td>
<td>
1805775
</td>
<td>
Number of aligned reads classified as mature miRNA strand
</td>
</tr>
<tr>
<td>
21
</td>
<td>
star miRNA
</td>
<td>
int
</td>
<td>
0
</td>
<td>
Number of aligned reads classified as star miRNA strand
</td>
</tr>
<tr>
<td>
22
</td>
<td>
precursor miRNA
</td>
<td>
int
</td>
<td>
762
</td>
<td>
Number of aligned reads classified as precursor miRNA, i.e. aligned to the region between the ends of miRNA primary transcript and mature (or star) strands.
</td>
</tr>
<tr>
<td>
23
</td>
<td>
miRNA loop
</td>
<td>
int
</td>
<td>
193
</td>
<td>
Number of aligned reads classified as miRNA loop, i.e., aligned to the 6 bases after the mature strand, where the loop between the mature strand and star strand would be.
</td>
</tr>
<tr>
<td>
24
</td>
<td>
unannotated miRNA
</td>
<td>
int
</td>
<td>
546
</td>
<td>
Number of aligned reads classified as unannotated miRNA feature
</td>
</tr>
<tr>
<td>
25
</td>
<td>
snoRNA
</td>
<td>
int
</td>
<td>
58823
</td>
<td>
Number of aligned reads classified as snoRNA, i.e. genes in UCSC knownGenes and kgXref where geneSymbol is in the form "SNOR*" or the description contains the string "small nucleolar RNA".
</td>
</tr>
<tr>
<td>
26
</td>
<td>
tRNA
</td>
<td>
int
</td>
<td>
789829
</td>
<td>
Number of aligned reads classified as tRNA, i.e. UCSC's RepeatMasker repeat class tRNA
</td>
</tr>
<tr>
<td>
27
</td>
<td>
rRNA
</td>
<td>
int
</td>
<td>
200471
</td>
<td>
Number of aligned reads classified as rRNA, i.e. UCSC's RepeatMasker repeat class rRNA
</td>
</tr>
<tr>
<td>
28
</td>
<td>
snRNA
</td>
<td>
int
</td>
<td>
16853
</td>
<td>
Number of aligned reads classified as snRNA, i.e. UCSC's RepeatMasker repeat class snRNA
</td>
</tr>
<tr>
<td>
39
</td>
<td>
scRNA
</td>
<td>
int
</td>
<td>
19532
</td>
<td>
Number of aligned reads classified as scRNA, i.e. UCSC's RepeatMasker repeat class scRNA
</td>
</tr>
<tr>
<td>
30
</td>
<td>
srpRNA
</td>
<td>
int
</td>
<td>
989
</td>
<td>
Number of aligned reads classified as srpRNA, i.e. UCSC's RepeatMasker repeat class srpRNA
</td>
</tr>
<tr>
<td>
31
</td>
<td>
Other RepeatMasker RNAs
</td>
<td>
int
</td>
<td>
67
</td>
<td>
Number of aligned reads classified as RNA, i.e. UCSC's RepeatMasker repeat class RNA
</td>
</tr>
<tr>
<td>
32
</td>
<td>
Non-Coding Exon
</td>
<td>
int
</td>
<td>
246408
</td>
<td>
Number of aligned reads classified as Non-coding Exon, possibly RNA, excluding snoRNAs reported above
</td>
</tr>
<tr>
<td>
33
</td>
<td>
3' UTR
</td>
<td>
int
</td>
<td>
20683
</td>
<td>
Number of aligned reads classified as 3' UTR
</td>
</tr>
<tr>
<td>
34
</td>
<td>
5' UTR
</td>
<td>
int
</td>
<td>
42743
</td>
<td>
Number of aligned reads classified as 5' UTR
</td>
</tr>
<tr>
<td>
35
</td>
<td>
Coding Exon
</td>
<td>
int
</td>
<td>
37942
</td>
<td>
Number of aligned reads classified as Coding Exon
</td>
</tr>
<tr>
<td>
36
</td>
<td>
Intron
</td>
<td>
int
</td>
<td>
76987
</td>
<td>
Number of aligned reads classified as Intron
</td>
</tr>
<tr>
<td>
37
</td>
<td>
LINE
</td>
<td>
int
</td>
<td>
7449
</td>
<td>
Number of aligned reads classified as LINE, i.e., UCSC RepeatMasker repeat class LINE
</td>
</tr>
<tr>
<td>
38
</td>
<td>
SINE
</td>
<td>
int
</td>
<td>
4610
</td>
<td>
Number of aligned reads classified as SINE, i.e., UCSC RepeatMasker repeat class SINE
</td>
</tr>
<tr>
<td>
39
</td>
<td>
LTR
</td>
<td>
int
</td>
<td>
6357
</td>
<td>
Number of aligned reads classified as LTR, i.e., UCSC RepeatMasker repeat class LTR
</td>
</tr>
<tr>
<td>
40
</td>
<td>
Satellite
</td>
<td>
int
</td>
<td>
325
</td>
<td>
Number of aligned reads classified as Satellite, i.e., UCSC RepeatMasker repeat class Satellite
</td>
</tr>
<tr>
<td>
41
</td>
<td>
RepeatMasker DNA
</td>
<td>
int
</td>
<td>
1320
</td>
<td>
Number of aligned reads classified as RepeatMasker DNA
</td>
</tr>
<tr>
<td>
42
</td>
<td>
RepeatMasker Low complexity
</td>
<td>
int
</td>
<td>
0
</td>
<td>
Number of aligned reads classified as RepeatMasker Low complexity
</td>
</tr>
<tr>
<td>
43
</td>
<td>
RepeatMasker Simple repeat
</td>
<td>
int
</td>
<td>
1620
</td>
<td>
Number of aligned reads classified as RepeatMasker Simple repeat
</td>
</tr>
<tr>
<td>
44
</td>
<td>
RepeatMasker Other
</td>
<td>
int
</td>
<td>
21
</td>
<td>
Number of aligned reads classified as RepeatMasker Other
</td>
</tr>
<tr>
<td>
45
</td>
<td>
RepeatMasker
Unknown
</td>
<td>
int
</td>
<td>
9
</td>
<td>
Number of aligned reads classified as RepeatMasker Unknown
</td>
</tr>
<tr>
<td>
46
</td>
<td>
Unknown
</td>
<td>
int
</td>
<td>
1537428
</td>
<td>
Number of aligned reads not classified
</td>
</tr>
<tr>
<td>
47
</td>
<td>
% Total miRNA
</td>
<td></td>
<td>
40.70%
</td>
<td>
Percentage of aligned reads mapped to miRNA species
</td>
</tr>
<tr>
<td>
48
</td>
<td>
% Crossmapped
miRNA
</td>
<td></td>
<td>
5.80%
</td>
<td>
Percentage of aligned reads crossmapped to miRNA species
</td>
</tr>
<tr>
<td>
49
</td>
<td>
% mature miRNA
</td>
<td></td>
<td>
34.87%
</td>
<td>
Percentage of aligned reads classified as mature miRNA strand
</td>
</tr>
<tr>
<td>
50
</td>
<td>
% star miRNA
</td>
<td></td>
<td>
0.00%
</td>
<td>
Percentage of aligned reads classified as star miRNA strand
</td>
</tr>
<tr>
<td>
51
</td>
<td>
% precursor miRNA
</td>
<td></td>
<td>
0.01%
</td>
<td>
Percentage of aligned reads classified as precursor miRNA
</td>
</tr>
<tr>
<td>
52
</td>
<td>
% miRNA loop
</td>
<td></td>
<td>
0.00%
</td>
<td>
Percentage of aligned reads classified as miRNA loop
</td>
</tr>
<tr>
<td>
53
</td>
<td>
% unannotated
miRNA
</td>
<td></td>
<td>
0.01%
</td>
<td>
Percentage of aligned reads classified as unannotated miRNA feature
</td>
</tr>
<tr>
<td>
54
</td>
<td>
% snoRNA
</td>
<td></td>
<td>
1.14%
</td>
<td>
Percentage of aligned reads classified as snoRNA
</td>
</tr>
<tr>
<td>
55
</td>
<td>
% tRNA
</td>
<td></td>
<td>
15.25%
</td>
<td>
Percentage of aligned reads classified as tRNA
</td>
</tr>
<tr>
<td>
56
</td>
<td>
% rRNA
</td>
<td></td>
<td>
3.87%
</td>
<td>
Percentage of aligned reads classified as rRNA
</td>
</tr>
<tr>
<td>
57
</td>
<td>
% snRNA
</td>
<td></td>
<td>
0.33%
</td>
<td>
Percentage of aligned reads classified as snRNA
</td>
</tr>
<tr>
<td>
58
</td>
<td>
% scRNA
</td>
<td></td>
<td>
0.38%
</td>
<td>
Percentage of aligned reads classified as scRNA
</td>
</tr>
<tr>
<td>
59
</td>
<td>
% srpRNA
</td>
<td></td>
<td>
0.02%
</td>
<td>
Percentage of aligned reads classified as srpRNA
</td>
</tr>
<tr>
<td>
60
</td>
<td>
% Other
RepeatMasker RNAs
</td>
<td></td>
<td>
0.00%
</td>
<td>
Percentage of aligned reads classified as RNA by RepeatMaskers
</td>
</tr>
<tr>
<td>
61
</td>
<td>
% Non-Coding Exon
</td>
<td></td>
<td>
4.76%
</td>
<td>
Percentage of aligned reads classified as Non-coding Exon, possibly RNA
</td>
</tr>
<tr>
<td>
62
</td>
<td>
% 3' UTR
</td>
<td></td>
<td>
0.40%
</td>
<td>
Percentage of aligned reads classified as 3' UTR
</td>
</tr>
<tr>
<td>
63
</td>
<td>
% 5' UTR
</td>
<td></td>
<td>
0.83%
</td>
<td>
Percentage of aligned reads classified as 5' UTR
</td>
</tr>
<tr>
<td>
64
</td>
<td>
% Coding Exon
</td>
<td></td>
<td>
0.73%
</td>
<td>
Percentage of aligned reads classified as Coding Exon
</td>
</tr>
<tr>
<td>
65
</td>
<td>
% Intron
</td>
<td></td>
<td>
1.49%
</td>
<td>
Percentage of aligned reads classified as Intron
</td>
</tr>
<tr>
<td>
66
</td>
<td>
% LINE
</td>
<td></td>
<td>
0.14%
</td>
<td>
Percentage of aligned reads classified as LINE
</td>
</tr>
<tr>
<td>
67
</td>
<td>
% SINE
</td>
<td></td>
<td>
0.09%
</td>
<td>
Percentage of aligned reads classified as SINE
</td>
</tr>
<tr>
<td>
68
</td>
<td>
% LTR
</td>
<td></td>
<td>
0.12%
</td>
<td>
Percentage of aligned reads classified as LTR
</td>
</tr>
<tr>
<td>
69
</td>
<td>
% Satellite
</td>
<td></td>
<td>
0.01%
</td>
<td>
Percentage of aligned reads classified as Satellite
</td>
</tr>
<tr>
<td>
70
</td>
<td>% RepeatMasker DNA</td>
<td></td>
<td>
0.03%
</td>
<td>
Percentage of aligned reads classified as RepeatMasker DNA
</td>
</tr>
<tr>
<td>
71
</td>
<td>
% RepeatMasker Low complexity
</td>
<td></td>
<td>
0.00%
</td>
<td>
Percentage of aligned reads classified as RepeatMasker Low complexity
</td>
</tr>
<tr>
<td>
72
</td>
<td>
% RepeatMasker
Simple repeat
</td>
<td></td>
<td>
0.03%
</td>
<td>
Percentage of aligned reads classified as RepeatMasker Simple repeat
</td>
</tr>
<tr>
<td>
73
</td>
<td>
% RepeatMasker
Other
</td>
<td></td>
<td>
0.00%
</td>
<td>
Percentage of aligned reads classified as RepeatMasker Other
</td>
</tr>
<tr>
<td>
74
</td>
<td>
% RepeatMasker
Unknown
</td>
<td></td>
<td>
0.00%
</td>
<td>
Percentage of aligned reads classified as RepeatMasker Unknown
</td>
</tr>
<tr>
<td>
75
</td>
<td>
% Unknown
</td>
<td></td>
<td>
29.69%
</td>
<td>
Percentage of aligned reads not classified
</td>
</tr>
</tbody>
</table>

## III. Feature Expression

 For each of the non-miRNA features summarized in alignment_stats.csv,
 files were created which lists the specific gene name expressed, and
 the expression level as a read count, sorted by expression.

 ### <feature_type>.txt

 where \<feature_type> is one of the features listed in
 alignment_stats.csv, e.g. rRNA

<table>
<thead>
<td><strong>Column</strong></td>
<td>
<strong>Title</strong>
</td>
<td><strong>Type</strong></td>
<td>
<strong>Example</strong>
</td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>Feature id</td>
<td>
str
</td>
<td>LSU-rRNA_Hsa</td>
<td>
id of feature
</td>
</tr>
<tr>
<td>
2
</td>
<td>read count</td>
<td>
int
</td>
<td>
149532
</td>
<td>Number of reads counted for this feature</td>
</tr>
</tbody>
</table>

## IV. miRNA Expression Matrix

 expn_matrix.txt

 miRNA by sample matrix of miRNA expression read counts

<table>
<thead>
<td><strong>Column</strong></td>
<td>
<strong>Title</strong>
</td>
<td><strong>Type</strong></td>
<td><strong>Example</strong></td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>
Gene
</td>
<td>
str
</td>
<td>hsa-let-7a-1</td>
<td>
Name of pre-miRNA
</td>
</tr>
<tr>
<td>
2
</td>
<td>sample</td>
<td>
int
</td>
<td>
131569
</td>
<td>Number of reads counted for the pre-miRNA species</td>
</tr>
</tbody>
</table>

 expn_matrix_norm.txt

 miRNA by sample matrix of miRNA expression normalized read counts

<table>
<thead>
<td><strong>Column</strong></td>
<td>
<strong>Title</strong>
</td>
<td><strong>Type</strong></td>
<td>
<strong>Example</strong>
</td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>
Gene
</td>
<td>
str
</td>
<td>
hsa-let-7a-1
</td>
<td>
Name of pre-miRNA
</td>
</tr>
<tr>
<td>
2
</td>
<td>sample</td>
<td>
float
</td>
<td>49106.433040</td>
<td>
Number of reads counted for the pre-miRNA species divided by total number of reads in expn_matrix_miRNA. txt
</td>
</tr>
</tbody>
</table>

 expn_matrix_norm_log.txt

 miRNA by sample matrix of miRNA expression normalized and logarized
 read counts

<table>
<thead>
<td><strong>Column</strong></td>
<td><strong>Title</strong></td>
<td><strong>Type</strong></td>
<td><strong>Example</strong></td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>
Gene
</td>
<td>
str
</td>
<td>hsa-let-7a-1</td>
<td>
Name of pre-miRNA
</td>
</tr>
<tr>
<td>
2
</td>
<td>sample</td>
<td>
float
</td>
<td>15.583624</td>
<td>
Log2 of the number of reads counted for the pre-miRNA species divided by total number of reads in expn_matrix_miRNA. txt. Note: Inf in the log file is represented by 0 as a convenience to allow for numerical processing by programs that don't recognize "Inf".
</td>
</tr>
</tbody>
</table>

 expn_matrix_mimat.txt

 miRNA by sample matrix of mature miRNA expression read counts

<table>
<thead>
<td><strong>Column</strong></td>
<td>
<strong>Title</strong>
</td>
<td><strong>Type</strong></td>
<td>
<strong>Example</strong>
</td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>
Gene
</td>
<td>
str
</td>
<td>
hsa-let-7a-2.
MIMAT0010195
</td>
<td>
mature miRNA MIMAT ID concatenated to the name of the appropriate pre-miRNA (with any final duplicate location modifier removed)
</td>
</tr>
<tr>
<td>
2
</td>
<td>
sample
</td>
<td>
int
</td>
<td>
5
</td>
<td>
Number of reads counted for the mature miRNA species
</td>
</tr>
</tbody>
</table>

 expn_matrix_mimat_norm.txt

 miRNA by sample matrix of mature miRNA expression normalized read
 counts

<table>
<thead>
<td><strong>Column</strong></td>
<td>
<strong>Title</strong>
</td>
<td><strong>Type</strong></td>
<td>
<strong>Example</strong>
</td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>
Gene
</td>
<td>
str
</td>
<td>
hsa-let-7a-2.
MIMAT0010195
</td>
<td>
mature miRNA MIMAT ID concatenated to the name of the appropriate pre-miRNA (with any final duplicate location modifier removed)
</td>
</tr>
<tr>
<td>
2
</td>
<td>sample
</td>
<td>
float
</td>
<td>
1.867497
</td>
<td>
Number of reads counted for the mature miRNA species divided by total number of reads in expn_matrix_mimat_miRNA.txt multiplied by 1,000,000 (reads per million - RPM)
</td>
</tr>
</tbody>
</table>

 expn_matrix_mimat_norm_log.txt

 miRNA by sample matrix of mature miRNA expression normalized and
 logarized read counts

<table>
<thead>
<td><strong>Column</strong></td>
<td><strong>Title</strong></td>
<td><strong>Type</strong></td>
<td><strong>Example</strong></td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>
Gene
</td>
<td>
str
</td>
<td>hsa-let-7a-2.MIMAT0010195
</td>
<td>
mature miRNA MIMAT ID concatenated to the name of the appropriate pre-miRNA (with any final duplicate location modifier removed)
</td>
</tr>
<tr>
<td>
2
</td>
<td>sample</td>
<td>
float
</td>
<td>
0.901106
</td>
<td>
Log2 of the number of reads counted for the mature miRNA species divided by total number of reads in expn_matrix_mimat_miRNA.txt multiplied by 1,000,000 (reads per million - RPM). Note: Inf in the log file is represented by 0 as a convenience to allow for numerical processing by programs that don't recognize "Inf".
</td>
</tr>
</tbody>
</table>

## V. QC

 filtered_taglengths.csv, chastity_taglengths.csv and
 softclip_taglengths.csv

<table>
<thead>
<td><strong>Column</strong></td>
<td>
<strong>Title</strong>
</td>
<td><strong>Type</strong></td>
<td>
<strong>Example</strong>
</td>
<td>
<strong>Description</strong>
</td>
</thead>
<tbody>
<tr>
<td>
1
</td>
<td>
taglen
</td>
<td>
int
</td>
<td>
15
</td>
<td>
tag length
</td>
</tr>
<tr>
<td>
2
</td>
<td>crossmapped</td>
<td>
float
</td>
<td>
0
</td>
<td>Percentage of aligned, chastity-failed or softclipped reads that are crossmapped to miRNA features</td>
</tr>
<tr>
<td>
3...29
</td>
<td>
feature type
</td>
<td>
float
</td>
<td>2.3870247476615933</td>
<td>
Percentage of aligned, chastity-failed or softclipped reads that are mapped to the feature
</td>
</tr>
</tbody>
</table>

 where feature_type is one of the features listed in
 alignment_stats.csv, e.g. rRNA

 ## VI. Other Output Files

 A bed directory is created for each sample. These directories contain
 text files in BED format showing coverage of all reads which pass the
 alignment filters. The text files can be used to generate wig coverage
 files. They are also sorted and filtered so that they can be used for
 peak discovery and novel miRNA prediction.

 Saturation.jpg graphs show the number of reads aligned to miRNAs, and
 the number of miRNAs covered by these reads; tags graphs the
 percentage of reads aligning to various types of non-coding RNAs, and
 the read lengths of these aligned reads; adapter graphs the read
 length distribution of all reads, including adapter-adapter dimers
 with no RNA.
