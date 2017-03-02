\name{ms_relation}
\alias{ms_relation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The relationship between mean and SD
}
\description{
This function is to show the relationship between mean and sD within biological replicates in each group.
}
\usage{
ms_relation(result, condition1, condition2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
A list giving the information of readscount and partial read alignment summary in each bam file
}
  \item{condition1}{
A vector that is a character string to give the name of one condition. eg:condition1="FTO".
}
  \item{condition2}{
A vector that is a character string to give the name of another condition. eg:condition2="WT".
}
}
\details{
This function plots the mean and sd relationship figure within IP or Input biological replicates. This function is supposed to be called through the \code{Trumpet_report}
  wrapper. It is better not to call this function directly.
}
\value{
The figure showing the relationship between mean and sd to identify abnormal biological replicates will be generated.
}
\examples{
# This function will be called by the main funciton Trumpet_report. So you just need to know how to use the main function Trumpet_report.

?Trumpet_report

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
