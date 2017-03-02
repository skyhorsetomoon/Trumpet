\name{read_distribute}
\alias{read_distribute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The reads count distribution in each gene
}
\description{
This function is to show the readscount distribution in each gene.
}
\usage{
read_distribute(result, condition1, condition2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
A list giving the information of readscount and partial read alignment summary in each bam file
}
  \item{condition1}{
A vector that is a character string to give the name of one condition. eg:condition1="FTO". Defult:condition1=NULL.
}
  \item{condition2}{
A vector that is a character string to give the name of another condition. eg:condition2="WT". Defult:condition2=NULL.
}
}
\details{
This function plots the reads coverage distribution of genomic features under different gene expression. This function is supposed to be called through the \code{Trumpet_report}
  wrapper. It is better not to call this function directly.
}
\value{
A figure showing the reads coverage distribution of genomic features under different gene expression will be generated.
}
\examples{
# This function will be called by the main funciton Trumpet_report. So you just need to know how to use the main function Trumpet_report.

?Trumpet_report

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
