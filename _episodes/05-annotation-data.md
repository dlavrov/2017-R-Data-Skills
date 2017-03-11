---
title: "Working with Annotation Data"
teaching: 60
exercises: 30
questions:
- "How can I use Genomic Ranges with real data?"
objectives:
- "Be able to retreive genomic regions "
- "Be able to group genomic regions"
- "Be able to append two data frames"
- "Be able to articulate what a `factor` is and how to convert between `factor` and `character`."
- "Be able to find basic properties of a data frames including size, class or type of the columns, names, and first few rows."
keypoints:
- "Use `cbind()` to add a new column to a data frame."
- "Use `rbind()` to add a new row to a data frame."
- "Remove rows from a data frame."
- "Use `na.omit()` to remove rows from a data frame with `NA` values."
- "Use `levels()` and `as.character()` to explore and manipulate factors"
- "Use `str()`, `nrow()`, `ncol()`, `dim()`, `colnames()`, `rownames()`, `head()` and `typeof()` to understand structure of the data frame"
- "Read in a csv file using `read.csv()`"
- "Understand `length()` of a data frame"
---



## GenomicFeatures and rtracklayer

Here we will learn about two Bioconductor packages for importing and working with external data. Both 
packages have different purposes and connect with GenomicRanges. The first, GenomicFeatures, is designed 
for working with transcript-based genomic annotations. The second, rtracklayer, is designed for importing 
and exporting annotation data into a variety of different formats. 

### GenomicFeatures and TranscriptDb objects

GenomicFeatures provides methods for creating and working with TranscriptDb objects. These TranscriptDb 
objects wrap annotation data in a way that allows genomic features, like genes, transcripts, exons, and 
coding sequences (CDS), to be extracted in a consistent way, regardless of the organism and origin of 
the annotation data.

> ## R Packages for Data
>
>  While it may sound strange to use an R package that contains data rather than R code, it’s actually 
> a clever and appropriate use of an R package. Bioconductor uses packages for many types of data, including 
> transcript and organism annotation data, experimental data, compressed reference genomes, and microarray 
> and SNP platform details. Packages are a terrific way to unite data from multiple sources into a single 
> easily loaded and explicitly versioned shared resource.
>
{: .callout}

Let’s start by installing GenomicFeatures and the transcript annotation package for mouse, _Mus musculus_ 
(We can check which annotation packages are available on the [Bioconductor annotation package page](http://www.bioconductor.org/packages/release/data/annotation/)):


~~~
library(BiocInstaller)
biocLite("GenomicFeatures")
biocLite("TxDb.Mmusculus.UCSC.mm10.ensGene")
~~~
{: .r}

Notice that all transcript annotation packages use the same naming scheme:
TxDb.<organism>.<annotation-source>.<annotation-version>. This annotation is for mouse genome version 
mm10 (Genome Reference Consortium version GRCm38), and the annotation comes from UCSC’s Ensembl track.

Now we can start working with the data:


~~~
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
txdb
~~~
{: .r}



~~~
TxDb object:
# Db type: TxDb
# Supporting package: GenomicFeatures
# Data source: UCSC
# Genome: mm10
# Organism: Mus musculus
# Taxonomy ID: 10090
# UCSC Table: ensGene
# UCSC Track: Ensembl Genes
# Resource URL: http://genome.ucsc.edu/
# Type of Gene ID: Ensembl gene ID
# Full dataset: yes
# miRBase build ID: NA
# transcript_nrow: 94647
# exon_nrow: 348801
# cds_nrow: 226312
# Db created by: GenomicFeatures package from Bioconductor
# Creation time: 2016-09-29 04:15:25 +0000 (Thu, 29 Sep 2016)
# GenomicFeatures version at creation time: 1.25.17
# RSQLite version at creation time: 1.0.0
# DBSCHEMAVERSION: 1.1
~~~
{: .output}



~~~
class(txdb)
~~~
{: .r}



~~~
[1] "TxDb"
attr(,"package")
[1] "GenomicFeatures"
~~~
{: .output}

Under the hood, the transcriptDb object represents a SQLite database contained inside this R package. We 
don’t need to know any SQLite to interact with and extract data from this object; the GenomicFeatures package 
provides all methods we’ll need.

#### Retrieving genomic regions using `genes()`, `transcripts()`, `exons()`, `cds()`, and `promoters()`

Suppose we wanted to access all gene regions in _Mus musculus_ (in this version of Ensembl annotation). 
There’s a simple accessor function for this, unsurprisingly named `genes()`:


~~~
mm_genes <- genes(txdb)
head(mm_genes)
~~~
{: .r}



~~~
GRanges object with 6 ranges and 1 metadata column:
                     seqnames                 ranges strand |
                        <Rle>              <IRanges>  <Rle> |
  ENSMUSG00000000001     chr3 [108107280, 108146146]      - |
  ENSMUSG00000000003     chrX [ 77837901,  77853623]      - |
  ENSMUSG00000000028    chr16 [ 18780447,  18811987]      - |
  ENSMUSG00000000031     chr7 [142575529, 142578143]      - |
  ENSMUSG00000000037     chrX [161117193, 161258213]      + |
  ENSMUSG00000000049    chr11 [108343354, 108414396]      + |
                                gene_id
                            <character>
  ENSMUSG00000000001 ENSMUSG00000000001
  ENSMUSG00000000003 ENSMUSG00000000003
  ENSMUSG00000000028 ENSMUSG00000000028
  ENSMUSG00000000031 ENSMUSG00000000031
  ENSMUSG00000000037 ENSMUSG00000000037
  ENSMUSG00000000049 ENSMUSG00000000049
  -------
  seqinfo: 66 sequences (1 circular) from mm10 genome
~~~
{: .output}



~~~
length(mm_genes)
~~~
{: .r}



~~~
[1] 39017
~~~
{: .output}

GenomicFeatures has other functions for retrieving all transcripts, exons, coding sequences (CDS), 
and promoters: `transcripts()`, `exons()`, `cds()`, and `promoters()`. Consult the documentation 
for this family of functions for extracting information from transcriptDb objects at `help(transcripts)`.

#### Grouping GRangesList objects by transcript, gene, etc.

It’s often more natural to work with a GRangesList object of these types of features grouped by some other type of feature than working with a massive GRanges list object of everything. For example, we might want to retrieve all exons grouped by transcript or gene:


~~~
mm_exons_by_tx <- exonsBy(txdb, by="tx")
mm_exons_by_gn <- exonsBy(txdb, by="gene")
length(mm_exons_by_tx)
~~~
{: .r}



~~~
[1] 94647
~~~
{: .output}



~~~
length(mm_exons_by_gn)
~~~
{: .r}



~~~
[1] 39017
~~~
{: .output}



~~~
head(mm_exons_by_gn)
~~~
{: .r}



~~~
GRangesList object of length 6:
$ENSMUSG00000000001 
GRanges object with 9 ranges and 2 metadata columns:
      seqnames                 ranges strand |   exon_id   exon_name
         <Rle>              <IRanges>  <Rle> | <integer> <character>
  [1]     chr3 [108107280, 108109316]      - |     67759        <NA>
  [2]     chr3 [108109403, 108109612]      - |     67760        <NA>
  [3]     chr3 [108111935, 108112088]      - |     67761        <NA>
  [4]     chr3 [108112473, 108112602]      - |     67762        <NA>
  [5]     chr3 [108115763, 108115891]      - |     67763        <NA>
  [6]     chr3 [108118301, 108118458]      - |     67764        <NA>
  [7]     chr3 [108123542, 108123683]      - |     67765        <NA>
  [8]     chr3 [108123795, 108123837]      - |     67766        <NA>
  [9]     chr3 [108145888, 108146146]      - |     67767        <NA>

$ENSMUSG00000000003 
GRanges object with 9 ranges and 2 metadata columns:
      seqnames               ranges strand | exon_id exon_name
  [1]     chrX [77837901, 77838114]      - |  343835      <NA>
  [2]     chrX [77837902, 77838114]      - |  343836      <NA>
  [3]     chrX [77841860, 77841911]      - |  343837      <NA>
  [4]     chrX [77842515, 77842616]      - |  343838      <NA>
  [5]     chrX [77842897, 77843007]      - |  343839      <NA>
  [6]     chrX [77845019, 77845086]      - |  343840      <NA>
  [7]     chrX [77847975, 77848114]      - |  343841      <NA>
  [8]     chrX [77853409, 77853530]      - |  343842      <NA>
  [9]     chrX [77853409, 77853623]      - |  343843      <NA>

...
<4 more elements>
-------
seqinfo: 66 sequences (1 circular) from mm10 genome
~~~
{: .output}

These functions that extract grouped features all take the transcriptDb object as their first argument 
and which type of feature to group by (e.g., gene, tx, exon, or cds) as their second argument. There 
are variety of these types of functions: `transcriptsBy()`, `exonsBy()`, `cdsBy()`, `intronsBy()`, 
`fiveUTRsByTranscript()`, and `threeUTRsByTranscript()` (see `help(transcriptsBy)` for more information).

#### Extracting subsets of features that overlap a specific chromosome

GenomicFeatures also provides functions for extracting subsets of features that overlap a specific 
chromosome or range. We can limit our queries to use a subset of chromosomes by setting which sequences 
our transcriptDb should query using the following approach:


~~~
seqlevels(txdb)
~~~
{: .r}



~~~
 [1] "chr1"                 "chr2"                 "chr3"                
 [4] "chr4"                 "chr5"                 "chr6"                
 [7] "chr7"                 "chr8"                 "chr9"                
[10] "chr10"                "chr11"                "chr12"               
[13] "chr13"                "chr14"                "chr15"               
[16] "chr16"                "chr17"                "chr18"               
[19] "chr19"                "chrX"                 "chrY"                
[22] "chrM"                 "chr1_GL456210_random" "chr1_GL456211_random"
[25] "chr1_GL456212_random" "chr1_GL456213_random" "chr1_GL456221_random"
[28] "chr4_GL456216_random" "chr4_GL456350_random" "chr4_JH584292_random"
[31] "chr4_JH584293_random" "chr4_JH584294_random" "chr4_JH584295_random"
[34] "chr5_GL456354_random" "chr5_JH584296_random" "chr5_JH584297_random"
[37] "chr5_JH584298_random" "chr5_JH584299_random" "chr7_GL456219_random"
[40] "chrX_GL456233_random" "chrY_JH584300_random" "chrY_JH584301_random"
[43] "chrY_JH584302_random" "chrY_JH584303_random" "chrUn_GL456239"      
[46] "chrUn_GL456359"       "chrUn_GL456360"       "chrUn_GL456366"      
[49] "chrUn_GL456367"       "chrUn_GL456368"       "chrUn_GL456370"      
[52] "chrUn_GL456372"       "chrUn_GL456378"       "chrUn_GL456379"      
[55] "chrUn_GL456381"       "chrUn_GL456382"       "chrUn_GL456383"      
[58] "chrUn_GL456385"       "chrUn_GL456387"       "chrUn_GL456389"      
[61] "chrUn_GL456390"       "chrUn_GL456392"       "chrUn_GL456393"      
[64] "chrUn_GL456394"       "chrUn_GL456396"       "chrUn_JH584304"      
~~~
{: .output}



~~~
seqlevels(txdb, force=TRUE) <- "chr1"
seqlevels(txdb)
~~~
{: .r}



~~~
[1] "chr1"
~~~
{: .output}



~~~
chr1_exons <- exonsBy(txdb, "tx")
all(unlist(seqnames(chr1_exons)) == "chr1")
~~~
{: .r}



~~~
Warning in setUnlistDataNames(x@unlistData, x@partitioning, use.names,
class(x)): failed to set names on the unlisted CompressedRleList object
~~~
{: .error}



~~~
[1] TRUE
~~~
{: .output}



~~~
txdb <- restoreSeqlevels(txdb) # restore txdb so it queries all sequences
~~~
{: .r}

#### Extracting subsets of features that overlap a specific chromosome

To extract feature data that only overlaps a specific region, use the following family of functions: 
`transcriptsByOverlaps()`, `exonsByOverlaps()`, and `cdsByOverlaps()` (see `help(transcriptByOverlaps()` 
for more information).

For example, say a QTL study has identified a quantitative trait loci in the region roughly on 
chromosome 8, from 123,260,562 to 123,557,264. Our coordinates are rough, so we’ll add 10kbp and 
get all genes within this expanded region:


~~~
qtl_region <- GRanges("chr8", IRanges(123260562, 123557264))
qtl_region
~~~
{: .r}



~~~
GRanges object with 1 range and 0 metadata columns:
      seqnames                 ranges strand
         <Rle>              <IRanges>  <Rle>
  [1]     chr8 [123260562, 123557264]      *
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
~~~
{: .output}



~~~
qtl_region_expanded <- qtl_region + 10e3
transcriptsByOverlaps(txdb, qtl_region_expanded)
~~~
{: .r}



~~~
GRanges object with 73 ranges and 2 metadata columns:
       seqnames                 ranges strand |     tx_id
          <Rle>              <IRanges>  <Rle> | <integer>
   [1]     chr8 [119910841, 124345722]      + |     47374
   [2]     chr8 [123254195, 123269745]      + |     47530
   [3]     chr8 [123254271, 123257636]      + |     47531
   [4]     chr8 [123254284, 123269743]      + |     47532
   [5]     chr8 [123254686, 123265070]      + |     47533
   ...      ...                    ...    ... .       ...
  [69]     chr8 [123559201, 123559319]      - |     49320
  [70]     chr8 [123560888, 123561006]      - |     49321
  [71]     chr8 [123562595, 123562713]      - |     49322
  [72]     chr8 [123564286, 123564404]      - |     49323
  [73]     chr8 [123565969, 123566087]      - |     49324
                  tx_name
              <character>
   [1] ENSMUST00000127664
   [2] ENSMUST00000001092
   [3] ENSMUST00000150356
   [4] ENSMUST00000156896
   [5] ENSMUST00000154450
   ...                ...
  [69] ENSMUST00000178208
  [70] ENSMUST00000179143
  [71] ENSMUST00000178297
  [72] ENSMUST00000179019
  [73] ENSMUST00000179081
  -------
  seqinfo: 66 sequences (1 circular) from mm10 genome
~~~
{: .output}

`transcriptByOverlaps()` returns all transcripts overlapping this range. All functions in this 
family also take a maxgap argument, which can be used to specify how large a gap between ranges 
is considered an overlap (0 by default). Setting the maxgap argument to 10kbp has the same effect 
as widening our ranges and then extracting elements as we did in the preceding example.

> ## Creating TranscriptDb Objects
>
> If there isn’t a transcript annotation package containing a transcriptDb object for your organism, 
> annotation track, or genome version of choice, GenomicFeatures provides a multitude of methods to 
> create one. If the annotation track you’d like to use is on the UCSC Genome Browser or a BioMart 
> (e.g., a data management system used by databases Ensembl and WormBase), GenomicFeatures contains the 
> functions make `TranscriptDbFromUCSC()` and `makeTranscriptDbFromBiomart()` for creating transcriptDb 
> from these databases. For some nonmodel systems, annotation only exists as a Gene Transfer Format (GTF) 
> or Gene Feature Format (GFF) file. In this case, a transcriptDb object can be created with the 
> `makeTranscriptDbFromGFF()` method. Once you’ve created a transcriptDb object using one of these methods, 
> you can save the underlying SQLite database with `saveDb()` and load it again with `loadDb()`. As a 
> demonstration of how easy this is, the following line downloads all required annotation data from Ensembl 
> to create a transcriptDb object for Platypus:
> ~~~
> species <- "oanatinus_gene_ensembl"
> platypus_txdb <- makeTranscriptDbFromBiomart("ensembl", species)
> ~~~
> Although creating transcriptDb directly and saving these objects as SQLite databases certainly works, 
> GenomicFeatures makes it easy to create a transcript annotation package directly from tracks from the 
> UCSC Genome Browser, Biomart, or from a transcriptDb object. See the functions `makeTxDbPackageFromUCSC()`, 
> `makeTxDbPackageFromBiomart()`, and `makeTxDbPackage()`, and the GenomicFeatures vignette or documentation 
> for more detail on how to use these functions.
>
{: .callout}

### rtracklayer package

#### Improting data into a GRanges object

The rtracklayer package includes flexible functions for importing and exporting data that stores 
ranges from a variety of formats like GTF/GFF, BED, BED Graph, and Wiggle. These functions automatically convert 
entries to GRanges objects and handle technicalities like missing values and assigning columns in the file to 
metadata columns—features that general solutions like read.delim() don’t have. Let’s look at how the rtracklayer 
function `import()` loads the Mus_musculus.GRCm38.75_chr1.gtf.gz file:


~~~
library(rtracklayer)
mm_gtf <- import('../data/Mus_musculus.GRCm38.75_chr1.gtf.gz')
colnames(mcols(mm_gtf)) # metadata columns read in
~~~
{: .r}



~~~
 [1] "source"            "type"              "score"            
 [4] "phase"             "gene_id"           "gene_name"        
 [7] "gene_source"       "gene_biotype"      "transcript_id"    
[10] "transcript_name"   "transcript_source" "tag"              
[13] "exon_number"       "exon_id"           "ccds_id"          
[16] "protein_id"       
~~~
{: .output}

The function `import()` detects the file type and imports all data as a GRanges object. There are also 
specific functions (e.g., `import.bed()`, `import.gff()`, `import.wig()`, etc.) that can you can use if you 
want to specify the format.

#### Exporting data

The rtracklayer package also provides export methods, for taking range data and saving it to a variety 
of common range formats. For example, suppose we wanted to write five random pseudogenes to a GTF file. 
We could use:


~~~
set.seed(0)
pseudogene_i <- which(mm_gtf$gene_biotype == "pseudogene" & mm_gtf$type == "gene")
pseudogene_sample <- sample(pseudogene_i, 5)
export(mm_gtf[pseudogene_sample], con="./export/five_random_pseudogene.gtf", format="GTF")
~~~
{: .r}

If we didn’t care about the specifics of these ranges (e.g., the information stored in the metadata 
columns), the BED file format may be more appropriate. BED files require at a minimum three columns: 
chromosomes (or sequence name), start position, and end position:


~~~
bed_data <- mm_gtf[pseudogene_sample]
mcols(bed_data) <- NULL # clear out metadata columns
export(bed_data, con="./export/five_random_pseudogene.bed", format="BED")
~~~
{: .r}


~~~
cat ./_episodes_rmd/export/five_random_pseudogene.bed
~~~
{: .r}




~~~
cat: ./_episodes_rmd/export/five_random_pseudogene.bed: No such file or directory
~~~
{: .output}

Finally, it’s worth noting that we’re just scratching the surface of rtracklayer’s capabilities. 
In addition to its import/export functions, rtracklayer also interfaces with genome browsers like 
UCSC’s Genome Browser. If you find yourself using the UCSC Genome Browser frequently, it’s worth 
reading the rtracklayer vignette and learning how to interact with it through R.

### Retrieving Promoter Regions: `flank()` and `promoters()`

Suppose we want to extract the promoter regions of all protein-coding genes from the GRCh38 _Mus musculus_ 
Ensembl GTF annotation track for chromosome 1 we loaded in using rtracklayer in the previous section 
(it contains additional information about the type of transcript (such as the gene_bio type and type columns). 
So, first we could find the subset of genes we’re interested in, let’s say all protein coding genes:


~~~
table(mm_gtf$gene_biotype) # Calling table() on gene_bioype column returns the number of features of each biotype.
~~~
{: .r}



~~~

             antisense                lincRNA                  miRNA 
                   480                    551                    354 
              misc_RNA polymorphic_pseudogene   processed_transcript 
                    93                     61                    400 
        protein_coding             pseudogene                   rRNA 
                 77603                    978                     69 
        sense_intronic      sense_overlapping                 snoRNA 
                    21                      4                    297 
                 snRNA 
                   315 
~~~
{: .output}



~~~
chr1_pcg <- mm_gtf[mm_gtf$type == "gene" & mm_gtf$gene_biotype == "protein_coding"] #subsetting 
summary(width(chr1_pcg))
~~~
{: .r}



~~~
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     78    9429   25750   60640   62420 1076000 
~~~
{: .output}



~~~
length(chr1_pcg)
~~~
{: .r}



~~~
[1] 1240
~~~
{: .output}



~~~
chr1_pcg_3kb_up <- flank(chr1_pcg, width=3000) #use flank to grab 3kbp upstream of each feature.
~~~
{: .r}

By default `flank()` takes strand into consideration (option ignore.strand=FALSE), so we just specify the 
width of our flanking region.

Extracting promoter regions is such a common operation that GenomicRanges packages have a convenience 
function to make it even simpler: `promoters()`. `promoters()` default arguments extract 3kbp upstream of 
each range, and 200bp downstream (but we can change this!):


~~~
chr1_pcg_3kb_up2 <- promoters(chr1_pcg, upstream=3000, downstream=0)
identical(chr1_pcg_3kb_up, chr1_pcg_3kb_up2)
~~~
{: .r}



~~~
[1] TRUE
~~~
{: .output}

## Retrieving Promoter Sequence: Connection GenomicRanges with Sequence Data

Once we’ve created promoter ranges using `flank()` or `promoters()`, we can use these to grab the promoter 
nucleotide sequences from a genome. We can do this entirely through Bioconductor's packages, but this approach 
requires that genome sequences be stored in a special R package. If your organism doesn’t have a premade genome 
package, it may be faster to write the promoter ranges to a file and use command-line BEDTools.


~~~
biocLite("BSgenome")
# Note: this file is about 712MB, so be ensure you have enough disk space before installing!
biocLite("BSgenome.Mmusculus.UCSC.mm10")
~~~
{: .r}

This is a BSgenome package, where BS stands for Biostrings, a Bioconductor package that contains classes 
for storing sequence data and methods for working with it. BSgenome packages contain the full reference 
genome for a particular organism, compressed and wrapped in a user-friendly package with common accessor 
methods. As always, it’s worth reading the vignettes for these packages on Bioconductor’s website.

Let’s first load the BSgenome.Mmusculus.UCSC.mm10 package and poke around:


~~~
library(BSgenome.Mmusculus.UCSC.mm10)
mm_gm <- BSgenome.Mmusculus.UCSC.mm10
organism(mm_gm)
~~~
{: .r}



~~~
[1] "Mus musculus"
~~~
{: .output}



~~~
providerVersion(mm_gm)
~~~
{: .r}



~~~
[1] "mm10"
~~~
{: .output}



~~~
provider(mm_gm)
~~~
{: .r}



~~~
[1] "UCSC"
~~~
{: .output}

We can use the accessor function `seqinfo()` to look at sequence information. BSgenome packages contain 
sequences for each chromosome, stored in a list-like structure we can access using indexing:


~~~
seqinfo(mm_gm)
~~~
{: .r}



~~~
Seqinfo object with 66 sequences (1 circular) from mm10 genome:
  seqnames       seqlengths isCircular genome
  chr1            195471971      FALSE   mm10
  chr2            182113224      FALSE   mm10
  chr3            160039680      FALSE   mm10
  chr4            156508116      FALSE   mm10
  chr5            151834684      FALSE   mm10
  ...                   ...        ...    ...
  chrUn_GL456392      23629      FALSE   mm10
  chrUn_GL456393      55711      FALSE   mm10
  chrUn_GL456394      24323      FALSE   mm10
  chrUn_GL456396      21240      FALSE   mm10
  chrUn_JH584304     114452      FALSE   mm10
~~~
{: .output}



~~~
mm_gm$chrM
~~~
{: .r}



~~~
  16299-letter "DNAString" instance
seq: GTTAATGTAGCTTAATAACAAAGCAAAGCACTGA...AATCATACTCTATTACGCAATAAACATTAACAA
~~~
{: .output}



~~~
mm_gm[[22]]
~~~
{: .r}



~~~
  16299-letter "DNAString" instance
seq: GTTAATGTAGCTTAATAACAAAGCAAAGCACTGA...AATCATACTCTATTACGCAATAAACATTAACAA
~~~
{: .output}

BSgenome objects can be searched using the string-matching and alignment functions in the Biostrings packages. 
These are meant for a few, quick queries (certainly not large-scale alignment or mapping!). For example:


~~~
library(Biostrings)
matchPattern("GGCGCGCC", mm_gm$chr1)
~~~
{: .r}



~~~
  Views on a 195471971-letter DNAString subject
subject: NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
views:
          start       end width
  [1]   4557138   4557145     8 [GGCGCGCC]
  [2]   4567326   4567333     8 [GGCGCGCC]
  [3]   6960128   6960135     8 [GGCGCGCC]
  [4]   7397441   7397448     8 [GGCGCGCC]
  [5]   7398352   7398359     8 [GGCGCGCC]
  ...       ...       ...   ... ...
[144] 191907520 191907527     8 [GGCGCGCC]
[145] 191934164 191934171     8 [GGCGCGCC]
[146] 191942448 191942455     8 [GGCGCGCC]
[147] 192834335 192834342     8 [GGCGCGCC]
[148] 193589224 193589231     8 [GGCGCGCC]
~~~
{: .output}

Using genomic sequence and _Mus musculus_ promoter regions we created in the previous section, 
we’re ready to extract promoter sequences. 

First, let's check that all sequences we want to grab are in the BSgenome object: 


~~~
all(seqlevels(chr1_pcg_3kb_up) %in% seqlevels(mm_gm))
~~~
{: .r}



~~~
[1] FALSE
~~~
{: .output}

Oops, because the BSgenome.Mmusculus.UCSC.mm10 file uses the UCSC chromosome name style and our 
annotation file uses Ensembl/NCBI style, this is not the case. Having to remap one chromosome naming 
scheme to another is quite a common operation as there is no standard chromosome naming scheme. 

First, let’s create a test GRanges object so we can show how we can manually change chromosome names. 
Bioconductor packages treat sequence names much like the factors. We can access and set the names of these 
levels using the seqlevels() function:


~~~
gr <- GRanges(c("chr1", "chr2"), IRanges(start=c(3, 4), width=10))
seqlevels(gr)
~~~
{: .r}



~~~
[1] "chr1" "chr2"
~~~
{: .output}



~~~
seqlevels(gr) <- c("1", "2")
seqlevels(gr)
~~~
{: .r}



~~~
[1] "1" "2"
~~~
{: .output}

Because having to switch between the style “chr1” (UCSC style) and “1” (Ensembl/NCBI style) is common, 
Bioconductor provides a convenience function `seqlevelsStyle()`:


~~~
seqlevelsStyle(chr1_pcg_3kb_up)
~~~
{: .r}



~~~
[1] "NCBI"    "Ensembl" "MSU6"    "AGPvF"  
~~~
{: .output}



~~~
seqlevelsStyle(mm_gm)
~~~
{: .r}



~~~
[1] "UCSC"
~~~
{: .output}



~~~
seqlevelsStyle(chr1_pcg_3kb_up) <- "UCSC"
all(seqlevels(chr1_pcg_3kb_up) %in% seqlevels(mm_gm))
~~~
{: .r}



~~~
[1] TRUE
~~~
{: .output}

With chromosome names consistent between our GRanges promoter regions and the mouse BSgenome package, 
it’s easy to grab the sequences for particular regions kept in a GRanges object:


~~~
chr1_3kb_seqs <- getSeq(mm_gm, chr1_pcg_3kb_up)
chr1_3kb_seqs
~~~
{: .r}



~~~
  A DNAStringSet instance of length 1240
       width seq
   [1]  3000 ATTCTGAGATGTGGTTACTAGATCAATGGG...CTAGCCGGGCCCAGCGCCCAGCCCCGCGG
   [2]  3000 GAAGTGGTATATCTGCCTAGTCTAGGTGTG...GTACTTAATCTGTGAGCACACATGCTAGT
   [3]  3000 CTTAAAAACCTAGATATTCTATTTTTTTTT...TGATAACGTCGTGAGCTCGGCTTCCAACA
   [4]  3000 GAATTGGCACAGTTTCACATGATTGGTCCA...CGGCCGCTGCAGCGCGACAGGGGCCGGGC
   [5]  3000 AAATATAAAGTTAACATACAAAAACTAGTC...GGGCGCGAGCTCGGGGCCGAACGCGAGGA
   ...   ... ...
[1236]  3000 CAACATGGGTAGTAGTGGGGGAGCTTTAGT...GGGCTGGCCTCACCAAGACGCAACAGGGA
[1237]  3000 AGGTGTGTTATATAATAATTGGTTTGACAC...AAAACTTGCTCTCTGGCTTCCTGGCGCCC
[1238]  3000 TTGGCCAGGTGATTGATCTTGTCCAACTGG...AGGCCGGGCTATATGCAAACCGAGTTCCC
[1239]  3000 GGCATTCCCCTATACTGGGGCATAGAACCT...TAAGGGTCTGCTCCCCACTGCTTACAGCC
[1240]  3000 GTAAATTTTCAGGTATATTTCTTTCTACTC...TGATATTTCTGTGGTCCTTATTTCTAGGT
~~~
{: .output}

We could then write these sequences to a FASTA file using:


~~~
writeXStringSet(chr1_3kb_seqs, file="./export/mm10_chr1_3kb_promoters.fasta", format="fasta")
~~~
{: .r}

It’s worth mentioning that Bioconductor has many other packages for working with promoter sequences, extracting motifs, and creating sequence logo plots. This functionality is well documented in a detailed workflow on the [Bioconductor website](http://www.bioconductor.org/help/workflows/gene-regulation-tfbs/).

## Getting Intergenic and Intronic Regions: `Gaps`, `Reduce`, and `Setdffs`

We already used `gaps()`.  Here we'll apply it to GRanges objects. It is important to consider
how range operations will work with strand and work across difference chromosomes/sequences. 
With IRanges, gaps were simple: they’re just the areas of a sequence with a range on them. With 
genomic ranges, gaps are calculated on every combination of strand and sequence:


~~~
gr2 <- GRanges(c("chr1", "chr2"), IRanges(start=c(4, 12), width=6), strand=c("+", "-"), seqlengths=c(chr1=21, chr2=41))
gr2 # so we can see what these ranges look like
~~~
{: .r}



~~~
GRanges object with 2 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1  [ 4,  9]      +
  [2]     chr2  [12, 17]      -
  -------
  seqinfo: 2 sequences from an unspecified genome
~~~
{: .output}


~~~
gaps(gr2)
~~~
{: .r}



~~~
GRanges object with 8 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1  [ 1,  3]      +
  [2]     chr1  [10, 21]      +
  [3]     chr1  [ 1, 21]      -
  [4]     chr1  [ 1, 21]      *
  [5]     chr2  [ 1, 41]      +
  [6]     chr2  [ 1, 11]      -
  [7]     chr2  [18, 41]      -
  [8]     chr2  [ 1, 41]      *
  -------
  seqinfo: 2 sequences from an unspecified genome
~~~
{: .output}

Were did all these gaps come from? When applied to GRanges, `gaps()` creates ranges for all sequences 
(chr1 and chr2 here) and all strands (+, -, and ambiguous strand *). For sequence-strand combinations 
without any ranges, the gap spans the entire chromosome.


~~~
gr3 <- gr2
strand(gr3) <- "*"
gr3
~~~
{: .r}



~~~
GRanges object with 2 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1  [ 4,  9]      *
  [2]     chr2  [12, 17]      *
  -------
  seqinfo: 2 sequences from an unspecified genome
~~~
{: .output}



~~~
gaps(gr3)[strand(gaps(gr3)) == "*"]
~~~
{: .r}



~~~
GRanges object with 4 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1  [ 1,  3]      *
  [2]     chr1  [10, 21]      *
  [3]     chr2  [ 1, 11]      *
  [4]     chr2  [18, 41]      *
  -------
  seqinfo: 2 sequences from an unspecified genome
~~~
{: .output}

Replacing strand with the ambiguous strand * is a common trick when we don’t care about keeping strand information.

### Extracting intergenic regions

Another approach to creating gaps using range operations is to use the set operations we discussed earlier. 
Let's use this approach to extract intergenic regions from all transcripts. This can be thought of as taking 
a set of ranges that represent entire chromosomes, and taking the set difference of these and all transcripts:


~~~
chrom_grngs <- as(seqinfo(txdb), "GRanges") #use the as() method to coerce the TranscriptDb object’s chromosome information into a GRanges object 
head(chrom_grngs, 2)
~~~
{: .r}



~~~
GRanges object with 2 ranges and 0 metadata columns:
       seqnames         ranges strand
          <Rle>      <IRanges>  <Rle>
  chr1     chr1 [1, 195471971]      *
  chr2     chr2 [1, 182113224]      *
  -------
  seqinfo: 66 sequences (1 circular) from mm10 genome
~~~
{: .output}



~~~
collapsed_tx <- reduce(transcripts(txdb)) #reduce transcripts from txdb, so overlapping transcripts are collapsed into a single range.
strand(collapsed_tx) <- "*"
intergenic <- setdiff(chrom_grngs, collapsed_tx) #take the set difference between ranges representing an entire chromosome, and those that represent transcripts on those ranges.
head(intergenic, 2)
~~~
{: .r}



~~~
GRanges object with 2 ranges and 0 metadata columns:
      seqnames             ranges strand
         <Rle>          <IRanges>  <Rle>
  [1]     chr1 [      1, 3054232]      *
  [2]     chr1 [3054734, 3102015]      *
  -------
  seqinfo: 66 sequences (1 circular) from mm10 genome
~~~
{: .output}

### Extracting introns

We'll now create GRanges objects representing the introns of transcripts. We’re going to do this two ways: 
first, using a simple convenience function appropriately named `intronsByTranscripts()`, then using range 
set operations.

#### Quick and easy approach

First, let’s consider the simple solution that uses the TranscriptDb object txdb we loaded earlier:


~~~
mm_introns <- intronsByTranscript(txdb)
head(mm_introns[['18880']], 2) # get first two introns for transcript 18880
~~~
{: .r}



~~~
GRanges object with 2 ranges and 0 metadata columns:
      seqnames                 ranges strand
         <Rle>              <IRanges>  <Rle>
  [1]     chr3 [113556174, 113558092]      -
  [2]     chr3 [113558219, 113558321]      -
  -------
  seqinfo: 66 sequences (1 circular) from mm10 genome
~~~
{: .output}

#### Manual approach

We’ll now look at the manual approach that teaches really important range manipulations. We’ll make this 
example simpler by only creating the introns for a single gene, amylase 1 (which from the Ensembl website 
has gene identifier ENSMUSG00000074264). The set operation approach we’ll take considers introns as the 
set difference between transcripts range and the exons’ ranges for these transcripts.


~~~
amy1 <- transcriptsBy(txdb, 'gene')$ENSMUSG00000074264
amy1
~~~
{: .r}



~~~
GRanges object with 5 ranges and 2 metadata columns:
      seqnames                 ranges strand |     tx_id
         <Rle>              <IRanges>  <Rle> | <integer>
  [1]     chr3 [113555710, 113577830]      - |     18879
  [2]     chr3 [113555953, 113574762]      - |     18880
  [3]     chr3 [113556149, 113562018]      - |     18881
  [4]     chr3 [113562690, 113574272]      - |     18882
  [5]     chr3 [113564987, 113606699]      - |     18883
                 tx_name
             <character>
  [1] ENSMUST00000067980
  [2] ENSMUST00000106540
  [3] ENSMUST00000172885
  [4] ENSMUST00000142505
  [5] ENSMUST00000174147
  -------
  seqinfo: 66 sequences (1 circular) from mm10 genome
~~~
{: .output}

Now, let's extract all exons from the TranscriptDb object. We will subset out the ones we need for these 
transcripts later:


~~~
mm_exons <- exonsBy(txdb, "tx")
mm_exons[[18881]]
~~~
{: .r}



~~~
GRanges object with 5 ranges and 3 metadata columns:
      seqnames                 ranges strand |   exon_id   exon_name
         <Rle>              <IRanges>  <Rle> | <integer> <character>
  [1]     chr3 [113561824, 113562018]      - |     68132        <NA>
  [2]     chr3 [113561632, 113561731]      - |     68130        <NA>
  [3]     chr3 [113558322, 113558440]      - |     68129        <NA>
  [4]     chr3 [113558093, 113558218]      - |     68128        <NA>
  [5]     chr3 [113556149, 113556173]      - |     68127        <NA>
      exon_rank
      <integer>
  [1]         1
  [2]         2
  [3]         3
  [4]         4
  [5]         5
  -------
  seqinfo: 66 sequences (1 circular) from mm10 genome
~~~
{: .output}

mm_exons contains all mouse exons, grouped by transcript, but our transcripts are in a single GRanges object 
(because these are transcripts grouped by gene, and we’re looking at one gene). So we split up transcript 
ranges into a GRangeList object, this time grouping by transcript identifier:


~~~
amy1_tx <- split(amy1, amy1$tx_id)
amy1_tx
~~~
{: .r}



~~~
GRangesList object of length 5:
$18879 
GRanges object with 1 range and 2 metadata columns:
      seqnames                 ranges strand |     tx_id
         <Rle>              <IRanges>  <Rle> | <integer>
  [1]     chr3 [113555710, 113577830]      - |     18879
                 tx_name
             <character>
  [1] ENSMUST00000067980

$18880 
GRanges object with 1 range and 2 metadata columns:
      seqnames                 ranges strand | tx_id            tx_name
  [1]     chr3 [113555953, 113574762]      - | 18880 ENSMUST00000106540

$18881 
GRanges object with 1 range and 2 metadata columns:
      seqnames                 ranges strand | tx_id            tx_name
  [1]     chr3 [113556149, 113562018]      - | 18881 ENSMUST00000172885

...
<2 more elements>
-------
seqinfo: 66 sequences (1 circular) from mm10 genome
~~~
{: .output}

Now we need to extract from mm_exons only the exons belonging to the amylase 1 gene 
transcripts. A nice feature of GRangesList objects created by the `split()` method is that each list element 
is given the name of the vector used to split the ranges. Similarly, each list element created by 
`exonsBy(txdb, "tx")` is also named, using the transcript names (because we used by="tx"). Thus, we can 
easily match up our transcripts and exons by transcript identifiers, which are the element names of both 
the GRangesList objects amy1_tx and mm_exons by using `match()`:


~~~
amy1_exons <- mm_exons[match(names(amy1_tx), names(mm_exons))]
amy1_exons
~~~
{: .r}



~~~
GRangesList object of length 5:
$18879 
GRanges object with 11 ranges and 3 metadata columns:
       seqnames                 ranges strand |   exon_id   exon_name
          <Rle>              <IRanges>  <Rle> | <integer> <character>
   [1]     chr3 [113577701, 113577830]      - |     68142        <NA>
   [2]     chr3 [113569845, 113570057]      - |     68139        <NA>
   [3]     chr3 [113569382, 113569528]      - |     68138        <NA>
   [4]     chr3 [113564869, 113565066]      - |     68136        <NA>
   [5]     chr3 [113563445, 113563675]      - |     68135        <NA>
   [6]     chr3 [113562628, 113562761]      - |     68133        <NA>
   [7]     chr3 [113561824, 113561946]      - |     68131        <NA>
   [8]     chr3 [113561632, 113561731]      - |     68130        <NA>
   [9]     chr3 [113558322, 113558440]      - |     68129        <NA>
  [10]     chr3 [113558093, 113558218]      - |     68128        <NA>
  [11]     chr3 [113555710, 113556173]      - |     68125        <NA>
       exon_rank
       <integer>
   [1]         1
   [2]         2
   [3]         3
   [4]         4
   [5]         5
   [6]         6
   [7]         7
   [8]         8
   [9]         9
  [10]        10
  [11]        11

...
<4 more elements>
-------
seqinfo: 66 sequences (1 circular) from mm10 genome
~~~
{: .output}

With both our exons and transcripts grouped by transcript, we can finally take the pairwise set 
difference (with `setdiff())` which creates the set of introns for each transcript. Remember, 
it’s imperative when using pairwise set functions to make sure your two objects are correctly matched up!


~~~
all(names(amy1_tx) == names(amy1_exons)) # check everything's matched up
~~~
{: .r}



~~~
[1] TRUE
~~~
{: .output}



~~~
amy1_introns <- setdiff(amy1_tx, amy1_exons) #psetdiff in the book did not work
head(amy1_introns[['18880']], 2) # the first two introns of amylase
~~~
{: .r}



~~~
GRanges object with 2 ranges and 0 metadata columns:
      seqnames                 ranges strand
         <Rle>              <IRanges>  <Rle>
  [1]     chr3 [113556174, 113558092]      -
  [2]     chr3 [113558219, 113558321]      -
  -------
  seqinfo: 66 sequences (1 circular) from mm10 genome
~~~
{: .output}

Are the introns we created manually for this gene identical to those created by the function intronsByTranscripts()?


~~~
identical(mm_introns[names(amy1_tx)], amy1_introns)
~~~
{: .r}



~~~
[1] TRUE
~~~
{: .output}

### Finding and Working with Overlapping Ranges

Finding and counting overlaps are probably the most important operations in working with GRanges 
objects and ranges in general. Bioconductor’s GenomicRanges package has functions for finding overlaps, 
subsetting by overlaps, and counting overlaps that are nearly identical to methods found in the IRanges 
package.

> ## Mind Your Overlaps (Part II)
>
> Overlaps get quite complex very quickly (as discussed in a warning earlier). Nowhere is this more 
> apparent than with RNA-seq, where many technical issues can make the simple act of estimating transcript 
> abundance (by counting how many aligned sequencing reads overlap a transcript region) incredibly 
> complicated. First, counting overlaps in RNA-seq yields quantitative results that are used in downstream 
> statistical analyses (e.g., the abundance of estimates for a particular transcript). This means that bias 
> and noise that enter the range overlap quantification process could lead to inaccuracies in differential 
> expression tests. Second, sequencing reads may align ambiguously—to multiple spots in the genome equally 
> well. Do we count these multiple mapping reads (sometimes known as “multireads”)? Or discard them entirely? 
> Some modern RNA-seq quantification methods like RSEM (Li et al., 2011) attempt to rescue these multireads. 
> Third, some reads may align uniquely, but overlap an exon that’s shared by two or more transcripts’ 
> isoforms. Should we ignore this read, count it once for each transcript (leading to double counting), 
> or assign it to a transcript? All of these rather technical decisions in finding and counting overlaps 
> can unfortunately lead to different biological results. Accurate RNA-seq quantification is still an actively
> researched area, and methods are still being developed and improving. In general, the methods in this 
> section are only appropriate for simple overlap operations, and may not be appropriate for quantification 
> tasks like estimating transcript abundance for RNA-seq.
>
{: .callout}

To demonstrate how `findOverlaps()` can be used with GRanges objects, we’ll load in a BED file of 
dbSNP (build 137) variants (in addition to SNPs, dbSNP also includes other types of variants like 
insertions/deletions, short tandem repeats, multinucleotide polymorphisms) for mouse chromosome 1 
(available in the Buffalo's book GitHub repository. Using rtracklayer, we’ll load these in:


~~~
library(rtracklayer)
dbsnp137 <- import("../data/mm10_snp137_chr1_trunc.bed.gz")
~~~
{: .r}

Suppose we want to find out how many of these variants fall into exonic regions, and how many do not. 
Using the mouse TranscriptDb object we loaded earlier (txdb) we can extract and collapse all overlapping 
exons with `reduce()`. We’ll also subset so that we’re only looking at chromosome 1 exons (because our 
variants are only from chromosome 1):


~~~
collapsed_exons <- reduce(exons(txdb), ignore.strand=TRUE)
chr1_collapsed_exons <- collapsed_exons[seqnames(collapsed_exons) == "chr1"]
~~~
{: .r}

Let’s explore our dbsnp137 object before looking for overlaps:


~~~
summary(width(dbsnp137))
~~~
{: .r}



~~~
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.000   1.000   1.000   1.142   1.000 732.000 
~~~
{: .output}



~~~
dbsnp137$name[which.max(width(dbsnp137))]
~~~
{: .r}



~~~
[1] "rs232497063"
~~~
{: .output}

The variant that’s 732 bases long is a bit large, so we find out which entry it is and grab its RS identifier. 
Checking on the [dbSNP website](https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=232497063), we see that 
this is indeed a real variant — an insertion/deletion on chromosome 1. From our `summary()`, we see 
that the vast majority of variants are 1 nucleotide long, which are either SNPs or 1 base pair insertons/deletions. 
Note that the minimum width is zero, hence, there are also variants with zero widths. Using dbsnp137[width(dbsnp137) == 0], 
we can take look at a few of these. In most cases, these correspond to insertions into the reference genome 
(these can be easily verified with the dbSNP or UCSC Genome Browser websites). Zero-width ranges will not overlap 
any feature, as they don’t have any region to overlap another range. To count these zero-width features too, we’ll 
resize using the `resize()` function:


~~~
dbsnp137_resized <- dbsnp137
zw_i <- width(dbsnp137_resized) == 0
dbsnp137_resized[zw_i] <- resize(dbsnp137_resized[zw_i], width=1)
~~~
{: .r}

With this set of ranges, it’s easy now to find out how many variants overlap our chromosome 1 exon regions. We’ll use findOverlaps() to create a Hits object. We’ll tell `findOverlaps()` to ignore strand:


~~~
hits <- findOverlaps(dbsnp137_resized, chr1_collapsed_exons, ignore.strand=TRUE)
hits
~~~
{: .r}



~~~
Hits object with 57624 hits and 0 metadata columns:
          queryHits subjectHits
          <integer>   <integer>
      [1]        89        2336
      [2]       190        1731
      [3]       266        9170
      [4]       356       11652
      [5]       426        5986
      ...       ...         ...
  [57620]   2699766       14422
  [57621]   2699875        9548
  [57622]   2699961        8735
  [57623]   2699985        7981
  [57624]   2699987        7691
  -------
  queryLength: 2700000 / subjectLength: 15048
~~~
{: .output}



~~~
length(unique(queryHits(hits)))
~~~
{: .r}



~~~
[1] 57623
~~~
{: .output}



~~~
length(unique(queryHits(hits)))/length(dbsnp137_resized)
~~~
{: .r}



~~~
[1] 0.02134185
~~~
{: .output}

Suppose we now wanted to look at the variants that do overlap these exons. We could do this by using the 
indices in the Hits object, but a simpler method is to use the method `subsetByOverlaps()`:


~~~
subsetByOverlaps(dbsnp137_resized, chr1_collapsed_exons, ignore.strand=TRUE)
~~~
{: .r}



~~~
GRanges object with 57623 ranges and 2 metadata columns:
          seqnames                 ranges strand |        name     score
             <Rle>              <IRanges>  <Rle> | <character> <numeric>
      [1]     chr1 [ 43032144,  43032144]      + | rs250123171         0
      [2]     chr1 [ 36713805,  36713805]      + |  rs50487270         0
      [3]     chr1 [132567494, 132567494]      + | rs247294715         0
      [4]     chr1 [160995431, 160995431]      + |  rs47617081         0
      [5]     chr1 [ 84036552,  84036553]      + | rs216202117         0
      ...      ...                    ...    ... .         ...       ...
  [57619]     chr1 [188263219, 188263219]      + |  rs13476293         0
  [57620]     chr1 [134780954, 134780954]      + | rs218301913         0
  [57621]     chr1 [130270464, 130270464]      + | rs266050681         0
  [57622]     chr1 [107380295, 107380295]      + | rs224267626         0
  [57623]     chr1 [ 98421207,  98421207]      + | rs224196900         0
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
~~~
{: .output}

Note that the length of this GRanges object matches up with the number of overlapping variants that overlap exons earlier.

GenomicRanges also includes a method for counting overlaps, `countOverlaps()`. So suppose we wanted to count 
the number of variants that overlap each exonic region. Because we want the counts to be based on the exonic 
regions (which were the subject ranges in these operations), we reverse the order of the arguments:


~~~
var_counts <- countOverlaps(chr1_collapsed_exons, dbsnp137_resized, ignore.strand=TRUE)
head(var_counts)
~~~
{: .r}



~~~
[1]  1  0 17 21  1  3
~~~
{: .output}

To make these easier to follow, let’s append them to our chromosome 1 exonic regions GRanges object as a metadata column:


~~~
chr1_collapsed_exons$num_vars <- var_counts
chr1_collapsed_exons
~~~
{: .r}



~~~
GRanges object with 15048 ranges and 1 metadata column:
          seqnames                 ranges strand |  num_vars
             <Rle>              <IRanges>  <Rle> | <integer>
      [1]     chr1     [3054233, 3054733]      * |         1
      [2]     chr1     [3102016, 3102125]      * |         0
      [3]     chr1     [3205901, 3207317]      * |        17
      [4]     chr1     [3213439, 3216968]      * |        21
      [5]     chr1     [3421702, 3421901]      * |         1
      ...      ...                    ...    ... .       ...
  [15044]     chr1 [195169702, 195169801]      * |         0
  [15045]     chr1 [195170991, 195171168]      * |         0
  [15046]     chr1 [195176553, 195176715]      * |         1
  [15047]     chr1 [195228278, 195228398]      * |         0
  [15048]     chr1 [195240910, 195241007]      * |         0
  -------
  seqinfo: 66 sequences (1 circular) from mm10 genome
~~~
{: .output}

You can check our analysis by  visiting some of these regions (like “chr1:3054233-3054733”) on the mouse UCSC Genome Browser

### Calculating Coverage of GRanges Objects

The same coverage methods we saw with IRanges also work with GRanges and GRangesList objects. Let’s 
generate some random fake 150bp reads on chromosome 19 (because it’s the smallest chromosome) of our 
mouse genome and then empirically measure their coverage. We’ll target 5x coverage. Using the famous 
Lander-Waterman coverage equation for coverage (C = LN/G, where C is coverage, L is the read length, 
N is the sequence length, and N is the number of reads), we see that for 150bp reads, a chromosome 
length of 61,431,566bp, and a target of 5x coverage, we need: 5*61,431,566/150 = 2,047,719 reads. 
Let’s generate these using R’s sample() function:


~~~
set.seed(0)
chr19_len <- seqlengths(txdb)['chr19']
chr19_len
~~~
{: .r}



~~~
   chr19 
61431566 
~~~
{: .output}



~~~
start_pos <- sample(1:(chr19_len-150), 2047719, replace=TRUE)
reads <- GRanges("chr19", IRanges(start=start_pos, width=150))
~~~
{: .r}

Now, let’s use the coverage() method from GenomicRanges to calculate the coverage of these random reads:


~~~
cov_reads <- coverage(reads)
cov_reads
~~~
{: .r}



~~~
RleList of length 1
$chr19
integer-Rle of length 61431559 with 3898230 runs
  Lengths: 61 17 35  4  7 46 33  8 10  7 ...  1  6 24 25  5 23 43 24 30 25
  Values :  0  1  2  3  4  5  6  7  6  7 ...  4  5  4  5  6  5  4  3  2  1
~~~
{: .output}

Coverage is calculated per every chromosome in the ranges object, and returned as a run-length encoded 
list (much like sequence names are returned by `seqnames()` from a GRangesList). We can calculate mean 
coverage per chromosome easily:


~~~
mean(cov_reads)
~~~
{: .r}



~~~
   chr19 
5.000001 
~~~
{: .output}

It’s also easy to calculate how much of this chromosome has no reads covering it (this will happen with 
shotgun sequencing, due to the random nature of read placement). We can do this two ways. First, we could 
use == and table():


~~~
table(cov_reads == 0)
~~~
{: .r}



~~~
         FALSE     TRUE
chr19 61025072   406487
~~~
{: .output}

Or, we could use some run-length encoding tricks (these are faster, and scale to larger data better):


~~~
sum(runLength(cov_reads)[runValue(cov_reads) == 0])
~~~
{: .r}



~~~
 chr19 
406487 
~~~
{: .output}



~~~
406487/chr19_len
~~~
{: .r}



~~~
      chr19 
0.006616908 
~~~
{: .output}

So, about 0.6% of our chromosome 19 remains uncovered. Interestingly, this is very close to the proportion of uncovered bases expected under Lander-Waterman, which is Poisson distributed with λ = c (where c is coverage).

> ## Challenge 3
>
> Read the output of `str(gapminder)` again;
> this time, use what you've learned about factors, lists and vectors,
> as well as the output of functions like `colnames` and `dim`
> to explain what everything that `str` prints out for gapminder means.
> If there are any parts you can't interpret, discuss with your neighbors!
>
> > ## Solution to Challenge 3
> >
> > The object `gapminder` is a data frame with columns
> > - `country` and `continent` are factors.
> > - `year` is an integer vector.
> > - `pop`, `lifeExp`, and `gdpPercap` are numeric vectors.
> >
> {: .solution}
{: .challenge}
