# Trumpet
**Transcriptome-guided quality assessment of m6A-seq data**
# About
Trumpet is designed for visualization of quality assessment of methylated RNA immunoprecipitation sequencing data. Trumpet package takes the aligned BAM files from m6A-seq data together with transcriptome information as the inputs to generate a quality assessment report in HTML format, which covers a number of metrics that are relevant to the m6A-seq technique. 
# Installation Instructions
User can dowload this package from Github by using the following method: Input the following commond in R or Rstudio
>install.packages("devtools")

>library(devtools)

> devtools::install_github("skyhorsetomoon/Trumpet")

> library(Trumpet)
# Usage Example
The following commands code will show how to use this package and generate the assessment report in HTML format:
**Case one**: Input the samples to be evaluated in BAM files and generate the HTML report directly.
### Input the bam file to assessment
> f1=system.file("extdata", "IP1.bam", package="Trumpet")

> f2=system.file("extdata", "IP2.bam", package="Trumpet")

> f3=system.file("extdata", "IP3.bam", package="Trumpet")

> f4=system.file("extdata", "IP4.bam", package="Trumpet")

> IP_BAM=c(f1,f2,f3,f4)

> f1=system.file("extdata", "Input1.bam", package="Trumpet")

> f2=system.file("extdata", "Input2.bam", package="Trumpet")

> f3=system.file("extdata", "Input3.bam", package="Trumpet")

> Input_BAM=c(f1,f2,f3)

> f1=system.file("extdata", "treated_IP1.bam", package="Trumpet")

> contrast_IP_BAM=c(f1)

> f2=system.file("extdata", "treated_Input1.bam", package="Trumpet")

> contrast_Input_BAM=c(f2)

### Input the annotation file

> GENE_ANNO_GTF <- system.file("extdata", "hg19toy.gtf", package="Trumpet")

### Generate the assessment report to call the main function **Trumpet\_report.R**

> trumpet_report <- Trumpet_report(IP_BAM,
                               Input_BAM,
                               contrast_IP_BAM,
                               contrast_Input_BAM,
                               condition1 = "untreated",
                               condition2 = "treated",
                               GENE_ANNO_GTF = GENE_ANNO_GTF)
                               
** Case two**: If user's Linux version can not generate HTML report directly, they can call the command **get\_readscount2.R** firstly to get the reads count saved as **.Rdata** format. Then, call the main function **Trumpet\_report.R** and set some parameters in Windows system. The following code show how to generate the report.
 ### Not to run the following commond
 > 
