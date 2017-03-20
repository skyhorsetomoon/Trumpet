\name{Trumpet_report}
\alias{Trumpet_report}
\title{
Transcriptome-guided quality assessment of methylated RNA immunoprecipitation sequencing data
}
\description{
The package is designed for visualization of quality assessment of Methylated RNA immunoprecipitation sequencing data.
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
               get_read=FALSE,
               INPUT_DIR=NA,
               OUTPUT_DIR=NA
	      )
}
\arguments{
  \item{IP_BAM}{
 a vector of file names or the concrete paths, which specifies a number of IP samples in bam format
}
  \item{Input_BAM}{
a vector of file names or the concrete paths, which specifies a number of Input samples in bam format
}
  \item{contrast_IP_BAM}{
a vector of file names or the concrete paths, which specifies a number of contrast IP samples in bam format
}
  \item{contrast_Input_BAM}{
a vector of file names or the concrete paths, which specifies a number of contrast Input samples in bam format
}
  \item{condition1}{
A vector that is a character string to give the name of one condition. eg:condition1="untreated".
}
  \item{condition2}{
A vector that is a character string to give the name of another condition. eg:condition2="treated".
}
  \item{GENE_ANNO_GTF}{
a string, which specifies a gene annotation GTF file if available, default: NA
}
  \item{TXDB}{
an optional TxDb object for gene annotation information used in the analysis, default: NA. Please refere to "GenomicFeatures" package for more details about the "TxDb" object.
}
  \item{sample_size}{
If you want the biological replicates readscount to be closer, you can set this parameter. eg:sample_size=10^8. Defult: sample_size=NA.
}
  \item{GENOME}{
a string,such as "hg19" or "mm10", which specifies the genome assembly used. If a gene annotation file is provided, this function will use it directly; otherwise, this function will 
 			download the gene annotation from UCSC using the genome assembly specified here and the gene annotation table specified in "UCSC_TABLE_NAME".
}
  \item{UCSC_TABLE_NAME}{
a string, which specifies the gene annotation used from UCSC, default: "knownGene". Please use function: supportedUCSCtables() to check available tables.
	Some tables may not be available for all genomes, and the "refGene" table does not work correctly due to multiple occuences of the same transcript on the same chromosome.
}
  \item{get_read}{
A logical variable, indicating whether to apply the function of get_readscount2 firstly. If it is TRUE, the commond of Trumpet_report will call the 'Trumpet_report2.Rmd' file and generate the quality assessment report in HTML format. If it is FALSE, this function will call the 'Trumpet_report.Rmd' file and generate the quality assessment report in HTML format. 
}
  \item{INPUT_DIR}{
A string, which specifies the directory of the "Input_Infor.Rdata" file, default: INPUT_DIR=NA. If user have the  "Input_Infor.Rdata" file getted from \code{get_readscount2} function, they shoul give the path of the file to INPUT_DIR.
After setting this parameter, user should \code{load(INPUT_DIR)} and then call the main function of \code{Trumpet_report}.
}
  \item{OUTPUT_DIR}{
A string, which specifies the output directory, default: OUTPUT_DIR=NA,the output result will save in the current directory. Otherwise, Trumpet will output the assessment result to the directory of user specified.
}
}
\details{
This function will output the reoprt of visualization of quality assessment of
    Methylated RNA immunoprecipitation sequencing data in HTML format.
}
\value{
This main function of the package will generate the MeRIP-seq quality assessment report in HTML format that can be saved in output directory setting by user.
}
\author{
Teng Zhang <2573552960@qq.com>
}
\examples{
# Input the bam file to assessment
f1=system.file("extdata", "IP1.bam", package="Trumpet")
f2=system.file("extdata", "IP2.bam", package="Trumpet")
f3=system.file("extdata", "IP3.bam", package="Trumpet")
f4=system.file("extdata", "IP4.bam", package="Trumpet")
IP_BAM=c(f1,f2,f3,f4)
f1=system.file("extdata", "Input1.bam", package="Trumpet")
f2=system.file("extdata", "Input2.bam", package="Trumpet")
f3=system.file("extdata", "Input3.bam", package="Trumpet")
Input_BAM=c(f1,f2,f3)
f1=system.file("extdata", "treated_IP1.bam", package="Trumpet")
contrast_IP_BAM=c(f1)
f2=system.file("extdata", "treated_Input1.bam", package="Trumpet")
contrast_Input_BAM=c(f2)
GENE_ANNO_GTF <- system.file("extdata", "hg19toy.gtf", package="Trumpet")
# Generate the assessment report
trumpet_report <- Trumpet_report(IP_BAM,
                               Input_BAM,
                               contrast_IP_BAM,
                               contrast_Input_BAM,
                               condition1 = "untreated",
                               condition2 = "treated",
                               GENE_ANNO_GTF = GENE_ANNO_GTF
			                          )
}
