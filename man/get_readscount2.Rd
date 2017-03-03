\name{get_readscount2}
\alias{get_readscount2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GetReadsCount2
}
\description{
This function is to get the mRNA readscount.
}
\usage{
get_readscount2(IP_BAM, 
		Input_BAM, 
		contrast_IP_BAM, 
		contrast_Input_BAM, 
		condition1, 
		condition2,
		GENE_ANNO_GTF = NA, 
		GENOME = NA,
		UCSC_TABLE_NAME = "knownGene",
		TXDB = NA, 
		sample_size = NA,  
		output_dir = NA
	       )
}
%- maybe also 'usage' for other objects documented here.
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
A vector that is a character string to give the name of another condition. eg:condition2="untreated".
}
  \item{condition2}{
A vector that is a character string to give the name of another condition. eg:condition2="treated".
}
  \item{GENE_ANNO_GTF}{
a string, which specifies a gene annotation GTF file if available, default: NA
}
  \item{GENOME}{
a string,such as "hg19" or "mm10", which specifies the genome assembly used. If a gene annotation file is provided, this function will use it directly; otherwise, this function will download the gene annotation from UCSC using the genome assembly specified here and the gene annotation table specified in "UCSC_TABLE_NAME".
}
  \item{UCSC_TABLE_NAME}{
a string, which specifies the gene annotation used from UCSC, default: "knownGene". Please use function: supportedUCSCtables() to check available tables.
	Some tables may not be available for all genomes, and the "refGene" table does not work correctly due to multiple occuences of the same transcript on the same chromosome.
}
  \item{TXDB}{
an optional TxDb object for gene annotation information used in the analysis, default: NA. Please refere to "GenomicFeatures" package for more details about the "TxDb" object.
}
  \item{sample_size}{
If you want the biological replicates readscount to be closer, you can set this parameter. eg:sample_size=10^8. Defult: sample_size=NA.
}
  \item{output_dir}{
a string, which specifies the output directory of the saving the reads count and other input parameter. If the output_dir=NA, the output result will save in the current directory. Otherwise, it will save the result to the directory of user specified.
}
}
\details{
This function will first get the reads count from BAM files and user can call the commond of Trumpet_report to get quality assessment report. If user's Liunx system can't generate HTML directly, they can call this commond firstly under Liunx system and the call the commond of Trumpet_report to generate HTML report under Windows system. 
}
\value{
This function will return a list that include two data frame having the read alignment summary, readscount in each bin and a table of input files information. The result will be saved in the "Input_Infor.Rdata" file in the path of output directory setting by user.
}	
\author{
Teng Zhang <2573552960@qq.com>
}
\examples{
#This function can help user get the readscount in each bin from the input BAM files and extract the read alignment summary. This function will save the result and some input paramter in Rdata format.
#This function is the basic module to assess Methylated RNA immunoprecipitation sequencing data.
# get two group BAM files and without sample eg:sample_size=NA by defult

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
# without sample method
result<-get_readscount2(IP_BAM,
                         Input_BAM,
                         contrast_IP_BAM,
                         contrast_Input_BAM,
                         condition1 = "untreated",
                         condition2 = "treated",
                         GENE_ANNO_GTF = GENE_ANNO_GTF 
		       )
# You can also use sample method that you can set sample_size=10^7.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
