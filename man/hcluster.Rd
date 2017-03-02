\name{hcluster}
\alias{hcluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The similarity between biological replicates
}
\description{
This function is to Show the similarity between  biological replicates in each group using hierarchical clustering method.
}
\usage{
hcluster(result)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
A list giving the information of readscount and partial read alignment summary in each bam file
}
}
\details{
In order to identify the similarity within the IP biological replicates or the Input biological replicates, we can use the hierachical cluster method and plot the result. This function is supposed to be called through the \code{Trumpet_report} wrapper. It is better not to call this function directly.
}
\value{
The plot will show the the similarity within the IP biological replicates or the Input biological replicates.
}
\examples{
# This function will be called by the main funciton Trumpet_report. So you just need to know how to use the main function Trumpet_report.

?Trumpet_report

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
