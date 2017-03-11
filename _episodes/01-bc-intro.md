---
title: "Introduction to Bioconductor"
teaching: 45
exercises: 10
questions:
- "How to find your way around RStudio?"
- "How to interact with R?"
- "How to manage your environment?"
- "How to install packages?"
objectives:
- "To gain familiarity with the various panes in the RStudio IDE"
- "To gain familiarity with the buttons, short cuts and options in the RStudio IDE"
- "To understand variables and how to assign to them"
- "To be able to manage your workspace in an interactive R session"
- "To be able to use mathematical and comparison operations"
- "To be able to call functions"
- "Introduction to package management"
keypoints:
- "Use RStudio to write and run R programs."
- "R has the usual arithmetic operators and mathematical functions."
- "Use `<-` to assign values to variables."
- "Use `ls()` to list the variables in a program."
- "Use `rm()` to delete objects in a program."
- "Use `install.packages()` to install packages (libraries)."
---



# Project Overview

## About

[Bioconductor](http://www.bioconductor.org/) is an open source and open development software project for the 
analysis of genome data (e.g. sequence, microarray, annotation and many other data types). 

- Statistical analysis: large data, technological artifacts, designed
  experiments; rigorous
- Comprehension: biological context, visualization, reproducibility
- High-throughput
    - Sequencing: RNASeq, ChIPSeq, variants, copy number, ...
    - Microarrays: expression, SNP, ...
    - Flow cytometry, proteomics, images, ...

## Finding help

The instructions for installing BioConductor packages are available in the administrative section of this manual. Documentation for Bioconductor packages can be found in the vignette of each package. A listing of the available packages is available on the Bioc Package page. Annotation libraries can be found here. The Bioc Mail GMANE Archive can be searched to find answers to questions that have been discussed on the mailing list. Another valuable information resource is the Bioconductor Book. The basic R help functions provide additional information about packages and their functions:

http://kasperdanielhansen.github.io/genbioconductor/
A lot of useful info!!!
git@github.com:kasperdanielhansen/genbioconductor.git

Packages, vignettes, work flows

- 1296 software packages; also...
    - 'Annotation' packages -- static data bases of identifier maps,
      gene models, pathways, etc; e.g., [TxDb.Hsapiens.UCSC.hg19.knownGene][]
    - 'Experiment packages -- data sets used to illustrate software
      functionality, e.g., [airway][]
- Discover and navigate via [biocViews][]
- Package 'landing page'
    - Title, author / maintainer, short description, citation,
      installation instructions, ..., download statistics
- All user-visible functions have help pages, most with runnable
  examples
- 'Vignettes' an important feature in Bioconductor -- narrative
  documents illustrating how to use the package, with integrated code
- 'Release' (every six months) and 'devel' branches
- [Support site](https://support.bioconductor.org);
  [videos](https://www.youtube.com/user/bioconductor), [recent
  courses](https://bioconductor.org/help/course-materials/)

Package installation and use

- A package needs to be installed once, using the instructions on the
  package landing page (e.g., [DESeq2][]).

    
    ```r
    source("https://bioconductor.org/biocLite.R")
    biocLite(c("DESeq2", "org.Hs.eg.db"))
    ```

- `biocLite()` installs _Bioconductor_, [CRAN][], and github packages.

- Once installed, the package can be loaded into an R session

    
    ```r
    library(GenomicRanges)
    ```

    and the help system queried interactively, as outlined above:

    
    ```r
    help(package="GenomicRanges")
    vignette(package="GenomicRanges")
    vignette(package="GenomicRanges", "GenomicRangesHOWTOs")
    ?GRanges
    ```

## Key concepts

Goals

- Reproducibility
- Interoperability
- Use

What a few lines of _R_ has to say


```r
x <- rnorm(1000)
y <- x + rnorm(1000)
df <- data.frame(X=x, Y=y)
plot(Y ~ X, df)
fit <- lm(Y ~ X, df)
anova(fit)
```

```
## Analysis of Variance Table
## 
## Response: Y
##            Df  Sum Sq Mean Sq F value    Pr(>F)    
## X           1 1098.99 1098.99  1118.7 < 2.2e-16 ***
## Residuals 998  980.43    0.98                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
abline(fit)
```

![plot of chunk five-lines](figure/five-lines-1.png)

Classes and methods -- "S3"

- `data.frame()`
  - Defines _class_ to coordinate data
  - Creates an _instance_ or _object_

- `plot()`, `lm()`, `anova()`, `abline()`: _methods_ defined on
  _generics_ to transform instances

- Discovery and help

    
    ```r
    class(fit)
    methods(class=class(fit))
    methods(plot)
    ?"plot"
    ?"plot.formula"
    ```

- tab completion!

_Bioconductor_ classes and methods -- "S4"

- Example: working with DNA sequences

    
    ```r
    library(Biostrings)
    dna <- DNAStringSet(c("AACAT", "GGCGCCT"))
    reverseComplement(dna)
    ```
    
    ```
    ##   A DNAStringSet instance of length 2
    ##     width seq
    ## [1]     5 ATGTT
    ## [2]     7 AGGCGCC
    ```

- Discovery and help

    
    ```r
    class(dna)
    ?"DNAStringSet-class"
    ?"reverseComplement,DNAStringSet-method"
    ```

## High-throughput sequence analysis work flows

1. Experimental design

2. Wet-lab sequence preparation (figure from http://rnaseq.uoregon.edu/)

    ![](our_figures/fig-rna-seq.png)

3. (Illumina) Sequencing (Bentley et al., 2008,
   doi:10.1038/nature07517)

    ![](http://www.nature.com/nature/journal/v456/n7218/images/nature07517-f1.2.jpg)

    - Primary output: FASTQ files of short reads and their [quality
      scores](http://en.wikipedia.org/wiki/FASTQ_format#Encoding)

4. Alignment
    - Choose to match task, e.g., [Rsubread][], Bowtie2 good for ChIPseq,
      some forms of RNAseq; BWA, GMAP better for variant calling
    - Primary output: BAM files of aligned reads
    - More recently: [kallisto][] and similar programs that produce
      tables of reads aligned to transcripts
5. Reduction
    - e.g., RNASeq 'count table' (simple spreadsheets), DNASeq called
      variants (VCF files), ChIPSeq peaks (BED, WIG files)
6. Analysis
    - Differential expression, peak identification, ...
7. Comprehension
    - Biological context

## _Bioconductor_ sequencing ecosystem

![Alt Sequencing Ecosystem](our_figures/SequencingEcosystem.png)

# High-Throughput Sequence Data

## DNA (and other) sequences

[Biostrings][]


```r
library(Biostrings)
data(phiX174Phage)
phiX174Phage
```

```
##   A DNAStringSet instance of length 6
##     width seq                                          names               
## [1]  5386 GAGTTTTATCGCTTCCATGAC...ATTGGCGTATCCAACCTGCA Genbank
## [2]  5386 GAGTTTTATCGCTTCCATGAC...ATTGGCGTATCCAACCTGCA RF70s
## [3]  5386 GAGTTTTATCGCTTCCATGAC...ATTGGCGTATCCAACCTGCA SS78
## [4]  5386 GAGTTTTATCGCTTCCATGAC...ATTGGCGTATCCAACCTGCA Bull
## [5]  5386 GAGTTTTATCGCTTCCATGAC...ATTGGCGTATCCAACCTGCA G97
## [6]  5386 GAGTTTTATCGCTTCCATGAC...ATTGGCGTATCCAACCTGCA NEB03
```

```r
letterFrequency(phiX174Phage, c("A", "C", "G", "T"))
```

```
##         A    C    G    T
## [1,] 1291 1157 1254 1684
## [2,] 1292 1156 1253 1685
## [3,] 1292 1156 1253 1685
## [4,] 1292 1155 1253 1686
## [5,] 1292 1156 1253 1685
## [6,] 1292 1155 1253 1686
```

```r
letterFrequency(phiX174Phage, "GC", as.prob=TRUE)
```

```
##            G|C
## [1,] 0.4476420
## [2,] 0.4472707
## [3,] 0.4472707
## [4,] 0.4470850
## [5,] 0.4472707
## [6,] 0.4470850
```

## Genomic ranges

[GenomicRanges][]

- `GRanges()`: genomic coordinates to represent annotations (exons,
  genes, regulatory marks, ...) and data (called peaks, variants,
  aligned reads)

    ![Alt GRanges](our_figures/GRanges.png)

- `GRangesList()`: genomic coordinates grouped into list elements
  (e.g., paired-end reads; exons grouped by transcript)

    ![Alt GRangesList](our_figures/GRangesList.png)

Operations

  ![Alt GRanges operations](our_figures/RangeOperations.png)

- intra-range: act on each range independently
    - e.g., `shift()`
- inter-range: act on all ranges in a `GRanges` object or
  `GRangesList` element
      - e.g., `reduce()`; `disjoin()`
- between-range: act on two separate `GRanges` or `GRangesList`
  objects
      - e.g., `findOverlaps()`, `nearest()`


```r
library(GenomicRanges)
gr <- GRanges("A", IRanges(c(10, 20, 22), width=5), "+")
shift(gr, 1)                            # intra-range
```

```
## GRanges object with 3 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]        A  [11, 15]      +
##   [2]        A  [21, 25]      +
##   [3]        A  [23, 27]      +
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
range(gr)                               # inter-range
```

```
## GRanges object with 1 range and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]        A  [10, 26]      +
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
reduce(gr)                              # inter-range
```

```
## GRanges object with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]        A  [10, 14]      +
##   [2]        A  [20, 26]      +
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
snps <- GRanges("A", IRanges(c(11, 17, 24), width=1))
findOverlaps(snps, gr)                  # between-range
```

```
## Hits object with 3 hits and 0 metadata columns:
##       queryHits subjectHits
##       <integer>   <integer>
##   [1]         1           1
##   [2]         3           2
##   [3]         3           3
##   -------
##   queryLength: 3 / subjectLength: 3
```

```r
setdiff(range(gr), gr)                  # 'introns'
```

```
## GRanges object with 1 range and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]        A  [15, 19]      +
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

## Summarized experiments

[SummarizedExperiment][]

![Alt SummarizedExperiment](our_figures/SummarizedExperiment.png)

- Coordinate feature x sample 'assays' with row (feature) and column
  (sample) descriptions.
- 'assays' can be any matrix-like object, including very large on-disk
  representations such as [HDF5Array][]

### _SummarizedExperiment_ exercise

The [airway][] experiment data package summarizes an RNA-seq
experiment investigating human smooth-muscle airway cell lines treated
with dexamethasone. Load the library and data set.


```r
library(airway)
data(airway)
airway
```

```
## class: RangedSummarizedExperiment 
## dim: 64102 8 
## metadata(1): ''
## assays(1): counts
## rownames(64102): ENSG00000000003 ENSG00000000005 ... LRG_98 LRG_99
## rowData names(0):
## colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
## colData names(9): SampleName cell ... Sample BioSample
```

`airway` is an example of the _SummarizedExperiment_ class. Explore
its `assay()` (the matrix of counts of reads overlapping genomic
regions of interest in each sample), `colData()` (a description of
each sample), and `rowRanges()` (a description of each region of
interest; here each region is an ENSEMBL gene).


```r
x <- assay(airway)
class(x)
```

```
## [1] "matrix"
```

```r
dim(x)
```

```
## [1] 64102     8
```

```r
head(x)
```

```
##                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
## ENSG00000000003        679        448        873        408       1138
## ENSG00000000005          0          0          0          0          0
## ENSG00000000419        467        515        621        365        587
## ENSG00000000457        260        211        263        164        245
## ENSG00000000460         60         55         40         35         78
## ENSG00000000938          0          0          2          0          1
##                 SRR1039517 SRR1039520 SRR1039521
## ENSG00000000003       1047        770        572
## ENSG00000000005          0          0          0
## ENSG00000000419        799        417        508
## ENSG00000000457        331        233        229
## ENSG00000000460         63         76         60
## ENSG00000000938          0          0          0
```

```r
colData(airway)
```

```
## DataFrame with 8 rows and 9 columns
##            SampleName     cell      dex    albut        Run avgLength
##              <factor> <factor> <factor> <factor>   <factor> <integer>
## SRR1039508 GSM1275862   N61311    untrt    untrt SRR1039508       126
## SRR1039509 GSM1275863   N61311      trt    untrt SRR1039509       126
## SRR1039512 GSM1275866  N052611    untrt    untrt SRR1039512       126
## SRR1039513 GSM1275867  N052611      trt    untrt SRR1039513        87
## SRR1039516 GSM1275870  N080611    untrt    untrt SRR1039516       120
## SRR1039517 GSM1275871  N080611      trt    untrt SRR1039517       126
## SRR1039520 GSM1275874  N061011    untrt    untrt SRR1039520       101
## SRR1039521 GSM1275875  N061011      trt    untrt SRR1039521        98
##            Experiment    Sample    BioSample
##              <factor>  <factor>     <factor>
## SRR1039508  SRX384345 SRS508568 SAMN02422669
## SRR1039509  SRX384346 SRS508567 SAMN02422675
## SRR1039512  SRX384349 SRS508571 SAMN02422678
## SRR1039513  SRX384350 SRS508572 SAMN02422670
## SRR1039516  SRX384353 SRS508575 SAMN02422682
## SRR1039517  SRX384354 SRS508576 SAMN02422673
## SRR1039520  SRX384357 SRS508579 SAMN02422683
## SRR1039521  SRX384358 SRS508580 SAMN02422677
```

```r
rowRanges(airway)
```

```
## GRangesList object of length 64102:
## $ENSG00000000003 
## GRanges object with 17 ranges and 2 metadata columns:
##        seqnames               ranges strand |   exon_id       exon_name
##           <Rle>            <IRanges>  <Rle> | <integer>     <character>
##    [1]        X [99883667, 99884983]      - |    667145 ENSE00001459322
##    [2]        X [99885756, 99885863]      - |    667146 ENSE00000868868
##    [3]        X [99887482, 99887565]      - |    667147 ENSE00000401072
##    [4]        X [99887538, 99887565]      - |    667148 ENSE00001849132
##    [5]        X [99888402, 99888536]      - |    667149 ENSE00003554016
##    ...      ...                  ...    ... .       ...             ...
##   [13]        X [99890555, 99890743]      - |    667156 ENSE00003512331
##   [14]        X [99891188, 99891686]      - |    667158 ENSE00001886883
##   [15]        X [99891605, 99891803]      - |    667159 ENSE00001855382
##   [16]        X [99891790, 99892101]      - |    667160 ENSE00001863395
##   [17]        X [99894942, 99894988]      - |    667161 ENSE00001828996
## 
## ...
## <64101 more elements>
## -------
## seqinfo: 722 sequences (1 circular) from an unspecified genome
```

It's easy to subset a _SummarizedExperiment_ on rows, columns and
assays, e.g., retaining just those samples in the `trt` level of the
`dex` factor. Accessing elements of the column data is common, so
there is a short-cut.


```r
cidx <- colData(airway)$dex %in% "trt"
airway[, cidx]
```

```
## class: RangedSummarizedExperiment 
## dim: 64102 4 
## metadata(1): ''
## assays(1): counts
## rownames(64102): ENSG00000000003 ENSG00000000005 ... LRG_98 LRG_99
## rowData names(0):
## colnames(4): SRR1039509 SRR1039513 SRR1039517 SRR1039521
## colData names(9): SampleName cell ... Sample BioSample
```

```r
## shortcut
airway[, airway$dex %in% "trt"]
```

```
## class: RangedSummarizedExperiment 
## dim: 64102 4 
## metadata(1): ''
## assays(1): counts
## rownames(64102): ENSG00000000003 ENSG00000000005 ... LRG_98 LRG_99
## rowData names(0):
## colnames(4): SRR1039509 SRR1039513 SRR1039517 SRR1039521
## colData names(9): SampleName cell ... Sample BioSample
```

It's also easy to perform range-based operations on
`SummarizedExperiment` objects, e.g., querying for range of chromosome
14 and then subsetting to contain only genes on this chromosome. Range
operations on rows are very common, so there are shortcuts here, too.


```r
chr14 <- as(seqinfo(rowRanges(airway)), "GRanges")["14"]
ridx <- rowRanges(airway) %over% chr14
airway[ridx,]
```

```
## class: RangedSummarizedExperiment 
## dim: 2244 8 
## metadata(1): ''
## assays(1): counts
## rownames(2244): ENSG00000006432 ENSG00000009830 ...
##   ENSG00000273259 ENSG00000273307
## rowData names(0):
## colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
## colData names(9): SampleName cell ... Sample BioSample
```

```r
## shortcut
chr14 <- as(seqinfo(airway), "GRanges")["14"]
airway[airway %over% chr14,]
```

```
## class: RangedSummarizedExperiment 
## dim: 2244 8 
## metadata(1): ''
## assays(1): counts
## rownames(2244): ENSG00000006432 ENSG00000009830 ...
##   ENSG00000273259 ENSG00000273307
## rowData names(0):
## colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
## colData names(9): SampleName cell ... Sample BioSample
```

Use the `assay()` and `rowSums()` function to remove all rows from the
`airway` object that have 0 reads overlapping all samples. Summarize
the library size (column sums of `assay()`) and plot a histogram of
the distribution of reads per feature of interest.

## BED, GFF, GTF, WIG import and export

Genome annotations: BED, WIG, GTF, etc. files. E.g., GTF:

- Component coordinates

        7   protein_coding  gene        27221129    27224842    .   -   . ...
        ...
        7   protein_coding  transcript  27221134    27224835    .   -   . ...
        7   protein_coding  exon        27224055    27224835    .   -   . ...
        7   protein_coding  CDS         27224055    27224763    .   -   0 ...
        7   protein_coding  start_codon 27224761    27224763    .   -   0 ...
        7   protein_coding  exon        27221134    27222647    .   -   . ...
        7   protein_coding  CDS         27222418    27222647    .   -   2 ...
        7   protein_coding  stop_codon  27222415    27222417    .   -   0 ...
        7   protein_coding  UTR         27224764    27224835    .   -   . ...
        7   protein_coding  UTR         27221134    27222414    .   -   . ...

- Annotations

        gene_id "ENSG00000005073"; gene_name "HOXA11"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
        ...
        ... transcript_id "ENST00000006015"; transcript_name "HOXA11-001"; transcript_source "ensembl_havana"; tag "CCDS"; ccds_id "CCDS5411";
        ... exon_number "1"; exon_id "ENSE00001147062";
        ... exon_number "1"; protein_id "ENSP00000006015";
        ... exon_number "1";
        ... exon_number "2"; exon_id "ENSE00002099557";
        ... exon_number "2"; protein_id "ENSP00000006015";
        ... exon_number "2";
        ...

[rtracklayer][]

- `import()`: import various formats to `GRanges` and similar instances
- `export()`: transform from `GRanges` and similar types to BED, GTF, ...
- Also, functions to interactively drive UCSC genome browser with data
  from _R_ / _Bioconductor_

## FASTQ files

Sequenced reads: FASTQ files

    @ERR127302.1703 HWI-EAS350_0441:1:1:1460:19184#0/1
    CCTGAGTGAAGCTGATCTTGATCTACGAAGAGAGATAGATCTTGATCGTCGAGGAGATGCTGACCTTGACCT
    +
    HHGHHGHHHHHHHHDGG<GDGGE@GDGGD<?B8??ADAD<BE@EE8EGDGA3CB85*,77@>>CE?=896=:
    @ERR127302.1704 HWI-EAS350_0441:1:1:1460:16861#0/1
    GCGGTATGCTGGAAGGTGCTCGAATGGAGAGCGCCAGCGCCCCGGCGCTGAGCCGCAGCCTCAGGTCCGCCC
    +
    DE?DD>ED4>EEE>DE8EEEDE8B?EB<@3;BA79?,881B?@73;1?########################

[ShortRead][]

- `readFastq()`: input
- `FastqStreamer()`: iterate through FASTQ files
- `FastqSampler()`: sample from FASTQ files, e.g., for quality assessment
- Functions for trimming and filters FASTQ files, QA assessment

## Aligned reads

Aligned reads: BAM files

- Header

        @HD     VN:1.0  SO:coordinate
        @SQ     SN:chr1 LN:249250621
        @SQ     SN:chr10        LN:135534747
        @SQ     SN:chr11        LN:135006516
        ...
        @SQ     SN:chrY LN:59373566
        @PG     ID:TopHat       VN:2.0.8b       CL:/home/hpages/tophat-2.0.8b.Linux_x86_64/tophat --mate-inner-dist 150 --solexa-quals --max-multihits 5 --no-discordant --no-mixed --coverage-search --microexon-search --library-type fr-unstranded --num-threads 2 --output-dir tophat2_out/ERR127306 /home/hpages/bowtie2-2.1.0/indexes/hg19 fastq/ERR127306_1.fastq fastq/ERR127306_2.fastq

- Alignments: ID, flag, alignment and mate

        ERR127306.7941162       403     chr14   19653689        3       72M             =       19652348        -1413  ...
        ERR127306.22648137      145     chr14   19653692        1       72M             =       19650044        -3720  ...
        ERR127306.933914        339     chr14   19653707        1       66M120N6M       =       19653686        -213   ...

- Alignments: sequence and quality

        ... GAATTGATCAGTCTCATCTGAGAGTAACTTTGTACCCATCACTGATTCCTTCTGAGACTGCCTCCACTTCCC        *'%%%%%#&&%''#'&%%%)&&%%$%%'%%'&*****$))$)'')'%)))&)%%%%$'%%%%&"))'')%))
        ... TTGATCAGTCTCATCTGAGAGTAACTTTGTACCCATCACTGATTCCTTCTGAGACTGCCTCCACTTCCCCAG        '**)****)*'*&*********('&)****&***(**')))())%)))&)))*')&***********)****
        ... TGAGAGTAACTTTGTACCCATCACTGATTCCTTCTGAGACTGCCTCCACTTCCCCAGCAGCCTCTGGTTTCT        '******&%)&)))&")')'')'*((******&)&'')'))$))'')&))$)**&&****************

- Alignments: Tags

        ... AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:72 YT:Z:UU NH:i:2  CC:Z:chr22      CP:i:16189276   HI:i:0
        ... AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:72 YT:Z:UU NH:i:3  CC:Z:=  CP:i:19921600   HI:i:0
        ... AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:4  MD:Z:72 YT:Z:UU XS:A:+  NH:i:3  CC:Z:=  CP:i:19921465   HI:i:0
        ... AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:4  MD:Z:72 YT:Z:UU XS:A:+  NH:i:2  CC:Z:chr22      CP:i:16189138   HI:i:0

[GenomicAligments][]

- `readGAlignments()`: Single-end reads
- `readGAlignmentPairs()`, `readGAlignmentsList()`: paired end reads

Working with large files

- `ScanBamParam()`: restrict input
- `BamFile(, yieldSize=)`: iteration
- [GenomicFiles][] provides useful helpers, e.g., `reduceByYield()`

### _GenomicAlignments_ exercise

The [RNAseqData.HNRNPC.bam.chr14][] package is an example of an
experiment data package. It contains a subset of BAM files used in a
gene knock-down experiment, as described in
`?RNAseqData.HNRNPC.bam.chr14`. Load the package and get the path to
the BAM files.


```r
library(RNAseqData.HNRNPC.bam.chr14)
fls = RNAseqData.HNRNPC.bam.chr14_BAMFILES
basename(fls)
```

```
## [1] "ERR127306_chr14.bam" "ERR127307_chr14.bam" "ERR127308_chr14.bam"
## [4] "ERR127309_chr14.bam" "ERR127302_chr14.bam" "ERR127303_chr14.bam"
## [7] "ERR127304_chr14.bam" "ERR127305_chr14.bam"
```

Create `BamFileList()`, basically telling R that these are paths to
BAM files rather than, say, text files from a spreadsheet.


```r
library(GenomicAlignments)
bfls = BamFileList(fls)
bfl = bfls[[1]]
```

Input and explore the aligments. See `?readGAlignments` and
`?GAlignments` for details on how to manipulate these objects.


```r
ga = readGAlignments(bfl)
ga
```

```
## GAlignments object with 800484 alignments and 0 metadata columns:
##            seqnames strand       cigar    qwidth     start       end
##               <Rle>  <Rle> <character> <integer> <integer> <integer>
##        [1]    chr14      +         72M        72  19069583  19069654
##        [2]    chr14      +         72M        72  19363738  19363809
##        [3]    chr14      -         72M        72  19363755  19363826
##        [4]    chr14      +         72M        72  19369799  19369870
##        [5]    chr14      -         72M        72  19369828  19369899
##        ...      ...    ...         ...       ...       ...       ...
##   [800480]    chr14      -         72M        72 106989780 106989851
##   [800481]    chr14      +         72M        72 106994763 106994834
##   [800482]    chr14      -         72M        72 106994819 106994890
##   [800483]    chr14      +         72M        72 107003080 107003151
##   [800484]    chr14      -         72M        72 107003171 107003242
##                width     njunc
##            <integer> <integer>
##        [1]        72         0
##        [2]        72         0
##        [3]        72         0
##        [4]        72         0
##        [5]        72         0
##        ...       ...       ...
##   [800480]        72         0
##   [800481]        72         0
##   [800482]        72         0
##   [800483]        72         0
##   [800484]        72         0
##   -------
##   seqinfo: 93 sequences from an unspecified genome
```

```r
table(strand(ga))
```

```
## 
##      +      -      * 
## 400242 400242      0
```

Many of the reads have cigar "72M". What does this mean? Can you
create a subset of reads that do not have this cigar? Interpret some
of the non-72M cigars. Any hint about what these cigars represent?


```r
tail(sort(table(cigar(ga))))
```

```
## 
## 18M123N54M 36M123N36M  64M316N8M 38M670N34M 35M123N37M        72M 
##        225        228        261        264        272     603939
```

```r
ga[cigar(ga) != "72M"]
```

```
## GAlignments object with 196545 alignments and 0 metadata columns:
##            seqnames strand       cigar    qwidth     start       end
##               <Rle>  <Rle> <character> <integer> <integer> <integer>
##        [1]    chr14      -     64M1I7M        72  19411677  19411747
##        [2]    chr14      + 55M2117N17M        72  19650072  19652260
##        [3]    chr14      - 43M2117N29M        72  19650084  19652272
##        [4]    chr14      - 40M2117N32M        72  19650087  19652275
##        [5]    chr14      + 38M2117N34M        72  19650089  19652277
##        ...      ...    ...         ...       ...       ...       ...
##   [196541]    chr14      -    51M1D21M        72 106950429 106950501
##   [196542]    chr14      +    31M1I40M        72 106960410 106960480
##   [196543]    chr14      +    52M1D20M        72 106965156 106965228
##   [196544]    chr14      -    13M1D59M        72 106965195 106965267
##   [196545]    chr14      -     6M1D66M        72 106965202 106965274
##                width     njunc
##            <integer> <integer>
##        [1]        71         0
##        [2]      2189         1
##        [3]      2189         1
##        [4]      2189         1
##        [5]      2189         1
##        ...       ...       ...
##   [196541]        73         0
##   [196542]        71         0
##   [196543]        73         0
##   [196544]        73         0
##   [196545]        73         0
##   -------
##   seqinfo: 93 sequences from an unspecified genome
```

Use the function `summarizeJunctions()` to identify genomic regions
that are spanned by reads with complicated cigars. Can you use the
argument `with.revmap=TRUE` to extract the reads supporting a
particular (e.g., first) junction?


```r
summarizeJunctions(ga)
```

```
## GRanges object with 4635 ranges and 3 metadata columns:
##          seqnames                 ranges strand |     score plus_score
##             <Rle>              <IRanges>  <Rle> | <integer>  <integer>
##      [1]    chr14   [19650127, 19652243]      * |         4          2
##      [2]    chr14   [19650127, 19653624]      * |         1          1
##      [3]    chr14   [19652355, 19653624]      * |         8          7
##      [4]    chr14   [19652355, 19653657]      * |         1          1
##      [5]    chr14   [19653773, 19653892]      * |         9          5
##      ...      ...                    ...    ... .       ...        ...
##   [4631]    chr14 [106912703, 106922227]      * |         1          0
##   [4632]    chr14 [106938165, 106938301]      * |        10          2
##   [4633]    chr14 [106938645, 106944774]      * |        24          7
##   [4634]    chr14 [106944969, 106950170]      * |         7          6
##   [4635]    chr14 [106950323, 106960260]      * |         1          1
##          minus_score
##            <integer>
##      [1]           2
##      [2]           0
##      [3]           1
##      [4]           0
##      [5]           4
##      ...         ...
##   [4631]           1
##   [4632]           8
##   [4633]          17
##   [4634]           1
##   [4635]           0
##   -------
##   seqinfo: 93 sequences from an unspecified genome
```

```r
junctions <- summarizeJunctions(ga, with.revmap=TRUE)
ga[ junctions$revmap[[1]] ]
```

```
## GAlignments object with 4 alignments and 0 metadata columns:
##       seqnames strand       cigar    qwidth     start       end     width
##          <Rle>  <Rle> <character> <integer> <integer> <integer> <integer>
##   [1]    chr14      + 55M2117N17M        72  19650072  19652260      2189
##   [2]    chr14      - 43M2117N29M        72  19650084  19652272      2189
##   [3]    chr14      - 40M2117N32M        72  19650087  19652275      2189
##   [4]    chr14      + 38M2117N34M        72  19650089  19652277      2189
##           njunc
##       <integer>
##   [1]         1
##   [2]         1
##   [3]         1
##   [4]         1
##   -------
##   seqinfo: 93 sequences from an unspecified genome
```

It is possible to do other actions on BAM files, e.g., calculating the
'coverage' (reads overlapping each base).


```r
coverage(bfl)$chr14
```

```
## integer-Rle of length 107349540 with 493510 runs
##   Lengths: 19069582       72   294083 ...       19       72   346298
##   Values :        0        1        0 ...        0        1        0
```
    
## Called variants: VCF files

- Header

          ##fileformat=VCFv4.2
          ##fileDate=20090805
          ##source=myImputationProgramV3.1
          ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
          ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
          ##phasing=partial
          ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
          ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
          ...
          ##FILTER=<ID=q10,Description="Quality below 10">
          ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
          ...
          ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
          ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">

- Location

          #CHROM POS     ID        REF    ALT     QUAL FILTER ...
          20     14370   rs6054257 G      A       29   PASS   ...
          20     17330   .         T      A       3    q10    ...
          20     1110696 rs6040355 A      G,T     67   PASS   ...

- Variant INFO

          #CHROM POS     ...	INFO                              ...
          20     14370   ...	NS=3;DP=14;AF=0.5;DB;H2           ...
          20     17330   ...	NS=3;DP=11;AF=0.017               ...
          20     1110696 ...	NS=2;DP=10;AF=0.333,0.667;AA=T;DB ...

- Genotype FORMAT and samples

          ... POS     ...  FORMAT      NA00001        NA00002        NA00003
          ... 14370   ...  GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
          ... 17330   ...  GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
          ... 1110696 ...  GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4

[VariantAnnotation][]

- `readVcf()`: VCF input
- `ScanVcfParam()`: restrict input to necessary fields / ranges
- `VcfFile()`: indexing and iterating through large VCF files
- `locateVariants()`: annotate in relation to genes, etc; see also
  [ensemblVEP][], [VariantFiltering][]
- `filterVcf()`: flexible filtering

[Bioconductor]: https://bioconductor.org
[CRAN]: https://cran.r-project.org
[biocViews]: https://bioconductor.org/packages/

[HDF5Array]: https://bioconductor.org/packages/HDF5Array
[AnnotationDbi]: https://bioconductor.org/packages/AnnotationDbi
[AnnotationHub]: https://bioconductor.org/packages/AnnotationHub
[BSgenome.Hsapiens.UCSC.hg19]: https://bioconductor.org/packages/BSgenome.Hsapiens.UCSC.hg19
[BSgenome]: https://bioconductor.org/packages/BSgenome
[BiocParallel]: https://bioconductor.org/packages/BiocParallel
[Biostrings]: https://bioconductor.org/packages/Biostrings
[CNTools]: https://bioconductor.org/packages/CNTools
[ChIPQC]: https://bioconductor.org/packages/ChIPQC
[ChIPseeker]: https://bioconductor.org/packages/ChIPseeker
[DESeq2]: https://bioconductor.org/packages/DESeq2
[DiffBind]: https://bioconductor.org/packages/DiffBind
[GenomicAlignments]: https://bioconductor.org/packages/GenomicAlignments
[GenomicFeatures]: https://bioconductor.org/packages/GenomicFeatures
[GenomicFiles]: https://bioconductor.org/packages/GenomicFiles
[GenomicRanges]: https://bioconductor.org/packages/GenomicRanges
[Gviz]: https://bioconductor.org/packages/Gviz
[Homo.sapiens]: https://bioconductor.org/packages/Homo.sapiens
[IRanges]: https://bioconductor.org/packages/IRanges
[KEGGREST]: https://bioconductor.org/packages/KEGGREST
[OmicCircos]: https://bioconductor.org/packages/OmicCircos
[PSICQUIC]: https://bioconductor.org/packages/PSICQUIC
[Rsamtools]: https://bioconductor.org/packages/Rsamtools
[Rsubread]: https://bioconductor.org/packages/Rsubread
[ShortRead]: https://bioconductor.org/packages/ShortRead
[SomaticSignatures]: https://bioconductor.org/packages/SomaticSignatures
[SummarizedExperiment]: https://bioconductor.org/packages/SummarizedExperiment
[TxDb.Hsapiens.UCSC.hg19.knownGene]: https://bioconductor.org/packages/TxDb.Hsapiens.UCSC.hg19.knownGene
[VariantAnnotation]: https://bioconductor.org/packages/VariantAnnotation
[VariantFiltering]: https://bioconductor.org/packages/VariantFiltering
[VariantTools]: https://bioconductor.org/packages/VariantTools
[airway]: https://bioconductor.org/packages/airway
[biomaRt]: https://bioconductor.org/packages/biomaRt
[cn.mops]: https://bioconductor.org/packages/cn.mops
[csaw]: https://bioconductor.org/packages/csaw
[edgeR]: https://bioconductor.org/packages/edgeR
[ensemblVEP]: https://bioconductor.org/packages/ensemblVEP
[epivizr]: https://bioconductor.org/packages/epivizr
[ggbio]: https://bioconductor.org/packages/ggbio
[h5vc]: https://bioconductor.org/packages/h5vc
[limma]: https://bioconductor.org/packages/limma
[metagenomeSeq]: https://bioconductor.org/packages/metagenomeSeq
[org.Hs.eg.db]: https://bioconductor.org/packages/org.Hs.eg.db
[org.Sc.sgd.db]: https://bioconductor.org/packages/org.Sc.sgd.db
[phyloseq]: https://bioconductor.org/packages/phyloseq
[rtracklayer]: https://bioconductor.org/packages/rtracklayer
[snpStats]: https://bioconductor.org/packages/snpStats

[dplyr]: https://cran.r-project.org/package=dplyr
[data.table]: https://cran.r-project.org/package=data.table
[Rcpp]: https://cran.r-project.org/package=Rcpp
[kallisto]: https://pachterlab.github.io/kallisto


> ## Challenge 1
>
> Which of the following are valid R variable names?
> 
> ```r
> min_height
> max.height
> _age
> .mass
> MaxLength
> min-length
> 2widths
> celsius2kelvin
> ```
>
> > ## Solution to challenge 1
> >
> > The following can be used as R variables:
> > 
> > ```r
> > min_height
> > max.height
> > MaxLength
> > celsius2kelvin
> > ```
> >
> > The following creates a hidden variable:
> > 
> > ```r
> > .mass
> > ```
> >
> > The following will not be able to be used to create a variable
> > 
> > ```r
> > _age
> > min-length
> > 2widths
> > ```
> {: .solution}
{: .challenge}

> ## Challenge 2
>
> What will be the value of each  variable  after each
> statement in the following program?
>
> 
> ```r
> mass <- 47.5
> age <- 122
> mass <- mass * 2.3
> age <- age - 20
> ```
>
> > ## Solution to challenge 2
> >
> > 
> > ```r
> > mass <- 47.5
> > ```
> > This will give a value of 47.5 for the variable mass
> >
> > 
> > ```r
> > age <- 122
> > ```
> > This will give a value of 122 for the variable age
> >
> > 
> > ```r
> > mass <- mass * 2.3
> > ```
> > This will multiply the existing value of 47.5 by 2.3 to give a new value of
> > 109.25 to the variable mass.
> >
> > 
> > ```r
> > age <- age - 20
> > ```
> > This will subtract 20 from the existing value of 122 to give a new value
> > of 102 to the variable age.
> {: .solution}
{: .challenge}


> ## Challenge 3
>
> Run the code from the previous challenge, and write a command to
> compare mass to age. Is mass larger than age?
>
> > ## Solution to challenge 3
> >
> > One way of answering this question in R is to use the `>` to set up the following:
> > 
> > ```r
> > mass > age
> > ```
> > 
> > ```
> > ## [1] TRUE
> > ```
> > This should yield a boolean value of TRUE since 109.25 is greater than 102.
> {: .solution}
{: .challenge}


> ## Challenge 4
>
> Clean up your working environment by deleting the mass and age
> variables.
>
> > ## Solution to challenge 4
> >
> > We can use the `rm` command to accomplish this task
> > 
> > ```r
> > rm(age, mass)
> > ```
> {: .solution}
{: .challenge}

> ## Challenge 5
>
> Install the following packages: `ggplot2`, `plyr`, `gapminder`
>
> > ## Solution to challenge 5
> >
> > We can use the `install.packages()` command to install the required packages.
> > 
> > ```r
> > install.packages("ggplot2")
> > install.packages("plyr")
> > install.packages("gapminder")
> > ```
> {: .solution}
{: .challenge}
