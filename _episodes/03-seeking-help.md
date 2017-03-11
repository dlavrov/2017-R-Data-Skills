---
title: "Seeking Help"
teaching: 10
exercises: 10
questions:
- "How can I get help in R?"
objectives:
- "To be able read R help files for functions and special operators."
- "To be able to use CRAN task views to identify packages to solve a problem."
- "To be able to seek help from your peers."
keypoints:
- "Use `help()` to get online help in R."
---



## Reading Help files

R, and every package, provide help files for functions. To search for help on a
function from a specific function that is in a package loaded into your
namespace (your interactive R session):


~~~
?function_name
help(function_name)
~~~
{: .r}

This will load up a help page in RStudio (or as plain text in R by itself).

Each help page is broken down into sections:

 - Description: An extended description of what the function does.
 - Usage: The arguments of the function and their default values.
 - Arguments: An explanation of the data each argument is expecting.
 - Details: Any important details to be aware of.
 - Value: The data the function returns.
 - See Also: Any related functions you might find useful.
 - Examples: Some examples for how to use the function.

Different functions might have different sections, but these are the main ones you should be aware of.

> ## Tip: Reading help files
>
> One of the most daunting aspects of R is the large number of functions
> available. It would be prohibitive, if not impossible to remember the
> correct usage for every function you use. Luckily, the help files
> mean you don't have to!
{: .callout}

## Special Operators

To seek help on special operators, use quotes:


~~~
?"+"
~~~
{: .r}

## Getting help on packages

Many packages come with "vignettes": tutorials and extended example documentation.
Without any arguments, `vignette()` will list all vignettes for all installed packages;
`vignette(package="package-name")` will list all available vignettes for
`package-name`, and `vignette("vignette-name")` will open the specified vignette.

If a package doesn't have any vignettes, you can usually find help by typing
`help("package-name")`.

## When you kind of remember the function

If you're not sure what package a function is in, or how it's specifically spelled you can do a fuzzy search:


~~~
??function_name
~~~
{: .r}

## When you have no idea where to begin

If you don't know what function or package you need to use
[CRAN Task Views](http://cran.at.r-project.org/web/views)
is a specially maintained list of packages grouped into
fields. This can be a good starting point.

## When your code doesn't work: seeking help from your peers

If you're having trouble using a function, 9 times out of 10,
the answers you are seeking have already been answered on
[Stack Overflow](http://stackoverflow.com/). You can search using
the `[r]` tag.

If you can't find the answer, there are a few useful functions to
help you ask a question from your peers:


~~~
?dput
~~~
{: .r}

Will dump the data you're working with into a format so that it can
be copy and pasted by anyone else into their R session.


~~~
sessionInfo()
~~~
{: .r}



~~~
R version 3.3.2 (2016-10-31)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X El Capitan 10.11.6

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    parallel  methods   stats     graphics  grDevices utils    
[8] datasets  base     

other attached packages:
 [1] ggplot2_2.2.1                      VariantAnnotation_1.20.2          
 [3] RNAseqData.HNRNPC.bam.chr14_0.12.0 ShortRead_1.32.0                  
 [5] GenomicAlignments_1.10.0           Rsamtools_1.26.1                  
 [7] BiocParallel_1.8.1                 rtracklayer_1.34.2                
 [9] airway_0.108.0                     SummarizedExperiment_1.4.0        
[11] Biobase_2.34.0                     GenomicRanges_1.26.2              
[13] GenomeInfoDb_1.10.3                Biostrings_2.42.1                 
[15] XVector_0.14.0                     IRanges_2.8.1                     
[17] S4Vectors_0.12.1                   BiocGenerics_0.20.0               
[19] checkpoint_0.3.18                  stringr_1.2.0                     
[21] knitr_1.15.1                      

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.9            plyr_1.8.4             highr_0.6             
 [4] RColorBrewer_1.1-2     GenomicFeatures_1.26.3 bitops_1.0-6          
 [7] tools_3.3.2            zlibbioc_1.20.0        biomaRt_2.30.0        
[10] digest_0.6.12          tibble_1.2             gtable_0.2.0          
[13] evaluate_0.10          RSQLite_1.1-2          memoise_1.0.0         
[16] lattice_0.20-34        BSgenome_1.42.0        Matrix_1.2-7.1        
[19] DBI_0.5-1              hwriter_1.3.2          grid_3.3.2            
[22] AnnotationDbi_1.36.2   XML_3.98-1.5           latticeExtra_0.6-28   
[25] magrittr_1.5           scales_0.4.1           assertthat_0.1        
[28] colorspace_1.3-2       labeling_0.3           stringi_1.1.2         
[31] lazyeval_0.2.0         munsell_0.4.3          RCurl_1.95-4.8        
~~~
{: .output}

Will print out your current version of R, as well as any packages you
have loaded. This can be useful for others to help reproduce and debug
your issue.

> ## Challenge 1
>
> Look at the help for the `c` function. What kind of vector do you
> expect you will create if you evaluate the following:
> 
> ~~~
> c(1, 2, 3)
> c('d', 'e', 'f')
> c(1, 2, 'f')
> ~~~
> {: .r}
> > ## Solution to Challenge 1
> >
> > The `c()` function creates a vector, in which all elements are the
> > same type. In the first case, the elements are numeric, in the
> > second, they are characters, and in the third they are characters:
> > the numeric values "coerced" to be characters.
> {: .solution}
{: .challenge}

> ## Challenge 2
>
> Look at the help for the `paste` function. You'll need to use this later.
> What is the difference between the `sep` and `collapse` arguments?
>
> > ## Solution to Challenge 2
> >
> > Look at the help for the `paste` function. You'll need to use this later.
> >
> > 
> > ~~~
> > help("paste")
> > ?paste
> > ~~~
> > {: .r}
> {: .solution}
{: .challenge}

> ## Challenge 3
> Use help to find a function (and its associated parameters) that you could
> use to load data from a csv file in which columns are delimited with "\t"
> (tab) and the decimal point is a "." (period). This check for decimal
> separator is important, especially if you are working with international
> colleagues, because different countries have different conventions for the
> decimal point (i.e. comma vs period).
> hint: use `??csv` to lookup csv related functions.
{: .challenge}

## Other ports of call

* [Quick R](http://www.statmethods.net/)
* [RStudio cheat sheets](http://www.rstudio.com/resources/cheatsheets/)
* [Cookbook for R](http://www.cookbook-r.com/)
