\name{get_readscount}
\alias{get_readscount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get reads count
}
\description{
This function is to get the mRNA readscount.
}
\usage{
get_readscount(IP_BAM, Input_BAM, contrast_IP_BAM, contrast_Input_BAM, GENE_ANNO_GTF = NA, GENOME = NA, UCSC_TABLE_NAME = "knownGene", TXDB = NA, sample_size = NA)
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
  \item{GENE_ANNO_GTF}{
 a string, which specifies a gene annotation GTF file if available, default: NA
}
  \item{GENOME}{
a string,such as "hg19" or "mm10", which specifies the genome assembly used. If a gene annotation file is provided, this function will use it directly; otherwise, this function will 
 			download the gene annotation from UCSC using the genome assembly specified here and the gene annotation table specified in "UCSC_TABLE_NAME".
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
}
\value{
This function will return a list that include two data frame having the read alignment summary, readscount in each bin and a table of input files information.
  \item{s }{a data frame including the readscount in each bin of every transcript}
  \item{read_alignment_summary}{a data frame including the read alignment summary}
  \item{transform_table}{a table transforming the input BAM files names or paths into another appropriate names}
}
\details{
This function is supposed to be called through the \code{Trumpet_report}
  wrapper. It is better not to call this function directly.
}
\examples{
# This function will be called by the main funciton Trumpet_report. So you just need to know how to use the main function Trumpet_report.

?Trumpet_report
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
