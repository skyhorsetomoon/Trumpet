\name{ctest_evluate}
\alias{ctest_evluate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assessment the enrichment signal using ctest method
}
\description{
This function is to evluate the enrichment signal in IP or reference_IP according to percentage of bins in different foldchange with ctest method

}
\usage{
ctest_evluate(result, condition1, condition2)
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
This function can calculate the percentages of enrichment signal using c-test method at different fold-change in IP samples compared with unified Input sample. This function is supposed to be called through the \code{Trumpet_report}
  wrapper. It is better not to call this function directly.
}
\value{
These figures show the calculation result about the percentages of enrichment signal using c-test method at different fold-change.
}
\examples{
# This function will be called by the main funciton Trumpet_report. So you just need to know how to use the main function Trumpet_report.
  
?Trumpet_report

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
