\name{Trumpet_report}
\alias{Trumpet_report}
\title{
Transcriptome-guided quality assessment of methylated RNA immunoprecipitation sequencing data
}
\description{
The package is designed for the quality assessment of methylated RNA immunoprecipitation sequencing data. 
}
\usage{
Trumpet_report(IP_BAM,
               Input_BAM,
               contrast_IP_BAM,
               contrast_Input_BAM,
               condition1,
               condition2,
               GENE_ANNO_GTF=NA,
               TXDB = NA,
               sample_size=NA,
               GENOME = NA,
               UCSC_TABLE_NAME = "knownGene",
               OUTPUT_DIR=NA)
}
\arguments{
  \item{IP_BAM}{
the path of all IP aligned MeRIP-seq data files in BAM format
}
  \item{Input_BAM}{
the path of all Input aligned MeRIP-seq data files in BAM format
}
  \item{contrast_IP_BAM}{
the path of all contrast IP aligned MeRIP-seq data files in BAM format
}
  \item{contrast_Input_BAM}{
the path of all contrast Input aligned MeRIP-seq data files in BAM format
}
  \item{condition1}{
A character string is used to give the name of one condition, e.g., condition1="untreated".
}
  \item{condition2}{
A character string is used to give the name of another condition. e.g., condition2="treated".
}
  \item{GENE_ANNO_GTF}{
A gene annotation file in GTF format if available
}
  \item{TXDB}{
An optional TxDb object for gene annotation information used in the analysis, default: NA. Please refere to "GenomicFeatures" package for more details about the "TxDb" object.
}
  \item{sample_size}{
To ensure the comparison is not affected by the difference in library size, user can set this parameter to ensure the same library size. eg:sample_size=10^8. Defult: sample_size=NA.
}
  \item{GENOME}{
A string,such as "hg19" or "mm10", which specifies the genome assembly used. If a gene annotation file is provided, this function will use it directly; otherwise, this function will download the gene annotation from UCSC using the genome assembly specified here and the gene annotation table specified in "UCSC_TABLE_NAME".
}
  \item{UCSC_TABLE_NAME}{
A string, which gives the gene annotation used from UCSC, default: "knownGene". Please use function: supportedUCSCtables() to check available tables.
	Some tables may not be available for all genomes, and the "refGene" table does not work correctly due to multiple occuences of the same transcript on the same chromosome.
}
  \item{OUTPUT_DIR}{
A string, which specify the output directory, default: OUTPUT_DIR=NA, the output result will save in the current directory. Otherwise, Trumpet will output the assessment report under the directory of user specified.
}
}
\details{
This function will output the quality assessment report of methylated RNA immunoprecipitation sequencing data in HTML format.
}
\value{
This main function of the package will generate the MeRIP-seq quality assessment report in HTML format that can be saved in output directory setting by user.
}
\author{
Teng Zhang <tengzhagn156@126.com>
}
\examples{
# Collect the path of all the aligned MeRIP-seq data files in BAM format.
f1 <- system.file("extdata", "IP1.bam", package="Trumpet")
f2 <- system.file("extdata", "IP2.bam", package="Trumpet")
f3 <- system.file("extdata", "IP3.bam", package="Trumpet")
f4 <- system.file("extdata", "IP4.bam", package="Trumpet")
f5 <- system.file("extdata", "Input1.bam", package="Trumpet")
f6 <- system.file("extdata", "Input2.bam", package="Trumpet")
f7 <- system.file("extdata", "Input3.bam", package="Trumpet")
f8 <- system.file("extdata", "treated_IP1.bam", package="Trumpet")
f9 <- system.file("extdata", "treated_Input1.bam", package="Trumpet")
ip_bam <- c(f1,f2,f3,f4)
input_bam <- c(f5,f6,f7)
contrast_ip_bam <- c(f8)
contrast_input_bam <- c(f9)
# Input the transcriptome annotation file
gtf <- system.file("extdata", "hg19toy.gtf", package="Trumpet")
# Generate the assessment report
# An HTML report will be generated under current working directory
trumpet_report <- Trumpet_report(IP_BAM = ip_bam,
                                   Input_BAM = input_bam,
                                   contrast_IP_BAM = contrast_ip_bam,
                                   contrast_Input_BAM = contrast_input_bam,
                                   condition1 = "untreated",
                                   condition2 = "treated",
                                   GENE_ANNO_GTF = gtf)
                               
# Browse the report                               
browseURL("Trumpet_report.html")                         
}
