---
title: "Genomic Ranges"
teaching: 40
exercises: 15
questions:
- "How can I read data in R?"
- "What are the basic data types in R?"
- "How do I represent categorical information in R?"
objectives:
- "To be aware of the different types of data."
- "To begin exploring data frames, and understand how it's related to vectors, factors and lists."
- "To be able to ask questions from R about the type, class, and structure of an object."
keypoints:
- "Use `read.csv` to read tabular data in R."
- "The basic data types in R are double, integer, complex, logical, and character."
- "Use factors to represent categories in R."
---



Because GenomicRanges extends IRanges, everything we’ve learned in the previous lesson can be directly 
applied to the genomic version of an IRanges object, GRanges. None of the function names nor behaviors 
differ much, besides two added complications: dealing with multiple chromosomes and strand. As we’ll 
see here, GenomicRanges manages these complications and greatly simplifies our lives when working with genomic data.

## Storing Genomic Ranges with GenomicRanges

The GenomicRanges package introduces a new class called GRanges for storing genomic ranges. The 
GRanges builds off of IRanges. IRanges objects are used to store ranges of genomic regions on a 
single sequence, and GRanges objects contain the two other pieces of information necessary to 
specify a genomic location: sequence name (e.g., which chromosome) and strand. GRanges objects 
also have metadata columns, which are the data linked to each genomic range. We can create GRanges 
objects much like we did with IRanges objects:


~~~
library(GenomicRanges)
gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"),
                    ranges=IRanges(start=5:8, width=10),
                    strand=c("+", "-", "-", "+"))
gr
~~~
{: .r}



~~~
GRanges object with 4 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1   [5, 14]      +
  [2]     chr1   [6, 15]      -
  [3]     chr2   [7, 16]      -
  [4]     chr3   [8, 17]      +
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
~~~
{: .output}

Using the GRanges() constructor, we can also add arbitrary metadata columns by specifying 
additional named arguments:


~~~
gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"), ranges=IRanges(start=5:8, width=10),
                    strand=c("+", "-", "-", "+"), gc=round(runif(4), 3))
gr
~~~
{: .r}



~~~
GRanges object with 4 ranges and 1 metadata column:
      seqnames    ranges strand |        gc
         <Rle> <IRanges>  <Rle> | <numeric>
  [1]     chr1   [5, 14]      + |     0.062
  [2]     chr1   [6, 15]      - |     0.206
  [3]     chr2   [7, 16]      - |     0.177
  [4]     chr3   [8, 17]      + |     0.687
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
~~~
{: .output}

This illustrates the structure of GRanges objects: genomic location specified by sequence name, 
range, and strand (on the left of the dividing bar), and metadata columns (on the right). 
Each row of metadata corresponds to a range on the same row.

All metadata attached to a GRanges object are stored in a DataFrame, which behaves identically to 
R’s base data.frame, but supports a wider variety of column types.

Whereas in the preceding example metadata columns are used to store numeric data, in practice we can 
store any type of data: identifiers and names (e.g., for genes, transcripts, SNPs, or exons), 
annotation data (e.g., conservation scores, GC content, repeat content, etc.), or experimental data 
(e.g., if ranges correspond to alignments, data like mapping quality and the number of gaps)

Also, notice seqlengths in the gr object we’ve just created. Because GRanges (and genomic range data 
in general) is always with respect to a particular genome version, we usually know beforehand what 
the length of each sequence/chromosome is.

We can specify the sequence lengths in the GRanges constructor, or set it after the object has been 
created using the seqlengths() function:


~~~
seqlens <- c(chr1=152, chr2=432, chr3=903)
gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"),
              ranges=IRanges(start=5:8, width=10),
              strand=c("+", "-", "-", "+"),
              gc=round(runif(4), 3),
              seqlengths=seqlens)
seqlengths(gr) <- seqlens # another way to do the same as above
gr
~~~
{: .r}



~~~
GRanges object with 4 ranges and 1 metadata column:
      seqnames    ranges strand |        gc
         <Rle> <IRanges>  <Rle> | <numeric>
  [1]     chr1   [5, 14]      + |     0.384
  [2]     chr1   [6, 15]      - |      0.77
  [3]     chr2   [7, 16]      - |     0.498
  [4]     chr3   [8, 17]      + |     0.718
  -------
  seqinfo: 3 sequences from an unspecified genome
~~~
{: .output}



~~~
seqlengths(gr)
~~~
{: .r}



~~~
chr1 chr2 chr3 
 152  432  903 
~~~
{: .output}

We access data in GRanges objects much like we access data from IRanges objects: with accessor functions:


~~~
# Same as with IRanges
start(gr)
~~~
{: .r}



~~~
[1] 5 6 7 8
~~~
{: .output}



~~~
end(gr)
~~~
{: .r}



~~~
[1] 14 15 16 17
~~~
{: .output}



~~~
width(gr)
~~~
{: .r}



~~~
[1] 10 10 10 10
~~~
{: .output}



~~~
# GRanges-specific
seqnames(gr)
~~~
{: .r}



~~~
factor-Rle of length 4 with 3 runs
  Lengths:    2    1    1
  Values : chr1 chr2 chr3
Levels(3): chr1 chr2 chr3
~~~
{: .output}



~~~
strand(gr)
~~~
{: .r}



~~~
factor-Rle of length 4 with 3 runs
  Lengths: 1 2 1
  Values : + - +
Levels(3): + - *
~~~
{: .output}

The returned objects are all run-length encoded. If we wish to extract all IRanges ranges from a GRanges 
object, we can use the ranges accessor function:


~~~
ranges(gr)
~~~
{: .r}



~~~
IRanges object with 4 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         5        14        10
  [2]         6        15        10
  [3]         7        16        10
  [4]         8        17        10
~~~
{: .output}

Like most objects in R, GRanges has a length that can be accessed with `length()`, and supports names:


~~~
length(gr)
~~~
{: .r}



~~~
[1] 4
~~~
{: .output}



~~~
names(gr) <- letters[1:length(gr)]
gr
~~~
{: .r}



~~~
GRanges object with 4 ranges and 1 metadata column:
    seqnames    ranges strand |        gc
       <Rle> <IRanges>  <Rle> | <numeric>
  a     chr1   [5, 14]      + |     0.384
  b     chr1   [6, 15]      - |      0.77
  c     chr2   [7, 16]      - |     0.498
  d     chr3   [8, 17]      + |     0.718
  -------
  seqinfo: 3 sequences from an unspecified genome
~~~
{: .output}

GRanges objects support the same style of subsetting as other R objects. For example, if you wanted all 
ranges with a start position greater than 7, you can use:


~~~
gr[start(gr) > 7] # where you use the T/F index created by `start(gr) > 7`:
~~~
{: .r}



~~~
GRanges object with 1 range and 1 metadata column:
    seqnames    ranges strand |        gc
       <Rle> <IRanges>  <Rle> | <numeric>
  d     chr3   [8, 17]      + |     0.718
  -------
  seqinfo: 3 sequences from an unspecified genome
~~~
{: .output}



~~~
start(gr) > 7
~~~
{: .r}



~~~
[1] FALSE FALSE FALSE  TRUE
~~~
{: .output}

Using the seqname() accessor, we can count how many ranges there are per chromosome and then subset to 
include only ranges for a particular chromosome:


~~~
table(seqnames(gr))
~~~
{: .r}



~~~

chr1 chr2 chr3 
   2    1    1 
~~~
{: .output}



~~~
gr[seqnames(gr) == "chr1"]
~~~
{: .r}



~~~
GRanges object with 2 ranges and 1 metadata column:
    seqnames    ranges strand |        gc
       <Rle> <IRanges>  <Rle> | <numeric>
  a     chr1   [5, 14]      + |     0.384
  b     chr1   [6, 15]      - |      0.77
  -------
  seqinfo: 3 sequences from an unspecified genome
~~~
{: .output}

The mcols() accessor is used access metadata columns:


~~~
mcols(gr)
~~~
{: .r}



~~~
DataFrame with 4 rows and 1 column
         gc
  <numeric>
1     0.384
2     0.770
3     0.498
4     0.718
~~~
{: .output}

Because this returns a DataFrame and DataFrame objects closely mimic data.frame, $ works to access 
specific columns. The usual syntactic shortcut for accessing a column works too:


~~~
mcols(gr)$gc
~~~
{: .r}



~~~
[1] 0.384 0.770 0.498 0.718
~~~
{: .output}



~~~
gr$gc
~~~
{: .r}



~~~
[1] 0.384 0.770 0.498 0.718
~~~
{: .output}

The real power of GRanges comes from combining subsetting with the data kept in our metadata columns. 
For example, we could easily compute the average GC content of all ranges on chr1:


~~~
mcols(gr[seqnames(gr) == "chr1"])$gc
~~~
{: .r}



~~~
[1] 0.384 0.770
~~~
{: .output}



~~~
mean(mcols(gr[seqnames(gr) == "chr1"])$gc)
~~~
{: .r}



~~~
[1] 0.577
~~~
{: .output}

If we wanted to find the average GC content for all chromosomes, we would use the same split-apply-combine 
strategy we learned before.

## Grouping Data with GRangesList

R’s lists can be used to group data together, such as after using split() to split a dataframe by a 
factor column. Grouping data this way is useful for both organizing data and processing it in chunks. 
GRanges objects also have their own version of a list, called GRangesList, which are similar to R’s 
lists. GRanges Lists can be created manually:


~~~
gr1 <- GRanges(c("chr1", "chr2"), IRanges(start=c(32, 95), width=c(24, 123)))
gr2 <- GRanges(c("chr8", "chr2"), IRanges(start=c(27, 12), width=c(42, 34)))
grl <- GRangesList(gr1, gr2)
grl
~~~
{: .r}



~~~
GRangesList object of length 2:
[[1]] 
GRanges object with 2 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1 [32,  55]      *
  [2]     chr2 [95, 217]      *

[[2]] 
GRanges object with 2 ranges and 0 metadata columns:
      seqnames   ranges strand
  [1]     chr8 [27, 68]      *
  [2]     chr2 [12, 45]      *

-------
seqinfo: 3 sequences from an unspecified genome; no seqlengths
~~~
{: .output}

GRangesList objects behave almost identically to R’s lists:


~~~
# unlist() combines all GRangesList elements into a single GRanges object:
unlist(grl)
~~~
{: .r}



~~~
GRanges object with 4 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1 [32,  55]      *
  [2]     chr2 [95, 217]      *
  [3]     chr8 [27,  68]      *
  [4]     chr2 [12,  45]      *
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
~~~
{: .output}



~~~
# we can combine many GRangesList objects with c():
doubled_grl <- c(grl, grl)
length(doubled_grl)
~~~
{: .r}



~~~
[1] 4
~~~
{: .output}

Like lists, we can also give and access list element names with the function names(). GRangesList 
objects also have some special features. For example, accessor functions for GRanges data (e.g., 
`seqnames()`, `start()`, `end()`, `width()`, `ranges()`, `strand()`, etc.) also work on GRangesList objects:


~~~
seqnames(grl)
~~~
{: .r}



~~~
RleList of length 2
[[1]]
factor-Rle of length 2 with 2 runs
  Lengths:    1    1
  Values : chr1 chr2
Levels(3): chr1 chr2 chr8

[[2]]
factor-Rle of length 2 with 2 runs
  Lengths:    1    1
  Values : chr8 chr2
Levels(3): chr1 chr2 chr8
~~~
{: .output}



~~~
start(grl)
~~~
{: .r}



~~~
IntegerList of length 2
[[1]] 32 95
[[2]] 27 12
~~~
{: .output}

Note the class of object Bioconductor uses for each of these: RleList and Integer List, which are
are analogous to GRangesList: a list for a specific type of data. RleList are lists for run-length 
encoded vectors, and IntegerList objects are lists for integers (with added features). 

In practice, we’re usually working with too much data to create GRanges objects manually with GRangesList(). 
More often, GRangesLists come about as the result of using the function split() on GRanges objects: 


~~~
chrs <- c("chr3", "chr1", "chr2", "chr2", "chr3", "chr1")
gr <- GRanges(chrs, IRanges(sample(1:100, 6, replace=TRUE),
                    width=sample(3:30, 6, replace=TRUE)))
gr
~~~
{: .r}



~~~
GRanges object with 6 ranges and 0 metadata columns:
      seqnames     ranges strand
         <Rle>  <IRanges>  <Rle>
  [1]     chr3 [100, 105]      *
  [2]     chr1 [ 39,  48]      *
  [3]     chr2 [ 78,  90]      *
  [4]     chr2 [ 94,  96]      *
  [5]     chr3 [ 22,  34]      *
  [6]     chr1 [ 66,  92]      *
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
~~~
{: .output}

And now splitting it:


~~~
gr_split <- split(gr, seqnames(gr))
gr_split[[1]]
~~~
{: .r}



~~~
GRanges object with 2 ranges and 0 metadata columns:
      seqnames     ranges strand
         <Rle>  <IRanges>  <Rle>
  [1]     chr3 [100, 105]      *
  [2]     chr3 [ 22,  34]      *
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
~~~
{: .output}



~~~
names(gr_split)
~~~
{: .r}



~~~
[1] "chr3" "chr1" "chr2"
~~~
{: .output}

Bioconductor also provides an `unsplit()` method to rejoin split data on the same factor that was 
used to split it. For example, because we created gr_split by splitting on seqnames(gr), we could 
unsplit gr_split with `unsplit(gr_split, seqnames(gr))`:


~~~
unsplit(gr_split, seqnames(gr))
~~~
{: .r}



~~~
GRanges object with 6 ranges and 0 metadata columns:
      seqnames     ranges strand
         <Rle>  <IRanges>  <Rle>
  [1]     chr3 [100, 105]      *
  [2]     chr1 [ 39,  48]      *
  [3]     chr2 [ 78,  90]      *
  [4]     chr2 [ 94,  96]      *
  [5]     chr3 [ 22,  34]      *
  [6]     chr1 [ 66,  92]      *
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
~~~
{: .output}

Grouped data is the basis of the split-apply-combine pattern. We can use `lapply()` and `sapply()` 
on GRangesLists objects to iterate through all elements and apply a function:


~~~
# Return the order of widths (smallest range to largest) of each GRanges element in a GRangesList.
lapply(gr_split, function(x) order(width(x)))
~~~
{: .r}



~~~
$chr3
[1] 1 2

$chr1
[1] 1 2

$chr2
[1] 2 1
~~~
{: .output}



~~~
# Return the start position of the earliest (leftmost) range:
sapply(gr_split, function(x) min(start(x)))
~~~
{: .r}



~~~
chr3 chr1 chr2 
  22   39   78 
~~~
{: .output}



~~~
# The number of ranges in every GRangesList object can be returned with this R idiom:
sapply(gr_split, length)
~~~
{: .r}



~~~
chr3 chr1 chr2 
   2    2    2 
~~~
{: .output}



~~~
# A faster approach to calculating element lengths is with the specialized function `elementLengths()`:
elementNROWS(gr_split)
~~~
{: .r}



~~~
chr3 chr1 chr2 
   2    2    2 
~~~
{: .output}

Although `lapply()` and `sapply()` give you the most freedom to write and use your own functions to 
apply to data. However, for many overlap operation functions (e.g., `reduce()`, `flank()`, `coverage()`, 
and `findOverlaps()`), we don’t need to explicitly apply them - they can work directly with GRangesList 
objects. For example, reduce() called on a GRangesList object automatically works at the list-element level:


~~~
reduce(gr_split)
~~~
{: .r}



~~~
GRangesList object of length 3:
$chr3 
GRanges object with 2 ranges and 0 metadata columns:
      seqnames     ranges strand
         <Rle>  <IRanges>  <Rle>
  [1]     chr3 [ 22,  34]      *
  [2]     chr3 [100, 105]      *

$chr1 
GRanges object with 2 ranges and 0 metadata columns:
      seqnames   ranges strand
  [1]     chr1 [39, 48]      *
  [2]     chr1 [66, 92]      *

$chr2 
GRanges object with 2 ranges and 0 metadata columns:
      seqnames   ranges strand
  [1]     chr2 [78, 90]      *
  [2]     chr2 [94, 96]      *

-------
seqinfo: 3 sequences from an unspecified genome; no seqlengths
~~~
{: .output}

reduce() illustrates an important (and extremely helpful) property of GRangesList objects: many methods 
applied to GRangesList objects work at the grouped-data level automatically.


> ## Challenge 7
>
> Consider the R output of the matrix below:
> 
> ~~~
>      [,1] [,2]
> [1,]    4    1
> [2,]    9    5
> [3,]   10    7
> ~~~
> {: .output}
> What was the correct command used to write this matrix? Examine
> each command and try to figure out the correct one before typing them.
> Think about what matrices the other commands will produce.
>
> 1. `matrix(c(4, 1, 9, 5, 10, 7), nrow = 3)`
> 2. `matrix(c(4, 9, 10, 1, 5, 7), ncol = 2, byrow = TRUE)`
> 3. `matrix(c(4, 9, 10, 1, 5, 7), nrow = 2)`
> 4. `matrix(c(4, 1, 9, 5, 10, 7), ncol = 2, byrow = TRUE)`
>
> > ## Solution to Challenge 7
> >
> > Consider the R output of the matrix below:
> > 
> > ~~~
> >      [,1] [,2]
> > [1,]    4    1
> > [2,]    9    5
> > [3,]   10    7
> > ~~~
> > {: .output}
> > What was the correct command used to write this matrix? Examine
> > each command and try to figure out the correct one before typing them.
> > Think about what matrices the other commands will produce.
> > 
> > ~~~
> > matrix(c(4, 1, 9, 5, 10, 7), ncol = 2, byrow = TRUE)
> > ~~~
> > {: .r}
> {: .solution}
{: .challenge}
