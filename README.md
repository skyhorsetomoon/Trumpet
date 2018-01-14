# Trumpet: an R package for transcriptome-guided quality assessment of m6A-seq data
Motivation: Methylated RNA immunoprecipitation sequencing (m6A-seq or MeRIP-seq) has been extensively used for profiling transcriptome-wide distribution of RNA N6-MethylAdenosine methyla-tion. However, due to the intrinsic properties of RNA molecule and the intricate procedures of this technique, m6A-seq data often suffers from various flaws. A convenient and comprehensive tool is solely needed to assess the quality of m6A-seq data to ensure it is suitable for subsequent analysis.

Results: From technical perspective, m6A-seq can be considered as a marriage of ChIP-seq and RNA-seq; hence, by effectively combing the data quality assessment metrics of the two techniques, we developed the Trumpet R package for evaluation of m6A-seq data quality. Trumpet package takes the aligned BAM files from m6A-seq data together with transcriptome information as the inputs to generate a quality assessment report in HTML format.

# Installation Instructions
The Trumpet package can be installed by the following R commands:
> devtools::install_github("skyhorsetomoon/Trumpet")
> 
> library(Trumpet)

# Usage Example
The following commands code will show how to use this package and generate the assessment report in HTML format.


## Input the samples to be evaluated in BAM files and generate the HTML report directly.

### Give all the bam files of MeRIP-Seq data
> f1 <- system.file("extdata", "IP1.bam", package="Trumpet")

> f2 <- system.file("extdata", "IP2.bam", package="Trumpet")

> f3 <- system.file("extdata", "IP3.bam", package="Trumpet")
> 
> f4 <- system.file("extdata", "IP4.bam", package="Trumpet")
> 
> f5 <- system.file("extdata", "Input1.bam", package="Trumpet")
> 
> f6 <- system.file("extdata", "Input2.bam", package="Trumpet")
> 
> f7 <- system.file("extdata", "Input3.bam", package="Trumpet")
> 
> f8 <- system.file("extdata", "treated_IP1.bam", package="Trumpet")
> 
> f9 <- system.file("extdata", "treated_Input1.bam", package="Trumpet")
> 
> ip\_bam <- c(f1,f2,f3,f4)
> 
> input\_bam <- c(f5,f6,f7)
> 
> contrast\_ip\_bam <- c(f8)
> 
> contrast\_input\_bam <- c(f9)

#### Input the annotation file
> gtf <- system.file("extdata", "hg19toy.gtf", package="Trumpet")

### Generate the assessment report 

> trumpet\_report <- Trumpet_report(IP_BAM = ip_bam,
>                                   Input\_BAM = input\_bam,
>                                   contrast\_IP_BAM = contrast\_ip\_bam,
>                                   contrast\_Input\_BAM = contrast\_input\_bam,
>                                   condition1 = "untreated",
>                                   condition2 = "treated",
>                                   sample\_size = 1*10^7,
>                                   GENE\_ANNO\_GTF = gtf)

### An HTML report will be generated under current working directory, which can be opened with the following command                        
> browseURL("Trumpet_report.html")
