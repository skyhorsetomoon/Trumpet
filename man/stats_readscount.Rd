\name{stats_readscount}
\alias{stats_readscount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The statistics of reads count percentage and read allignment summary
}
\description{
This function is to get the read allignment summary and the statistics of readscount percentage in each transcript and each bin.
}
\usage{
stats_readscount(result)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
A list giving the information of readscount and partial read alignment summary in each bam file 
}
}
\value{
  \item{transform_table }{A table transforming the input sample path into name }
  \item{readscount_summary }{A data frame including readscount summary information }
  \item{p_geom}{A dataframe including some percentages of readscount in each transcript at IP bam files}
  \item{pt_geom}{A dataframe including some percentages of readscount in each transcript at Input bam files}
  \item{p_bin}{A dataframe including some percentages of readscount in each transcript at IP bam files}
  \item{pt_bin}{A dataframe including some percentages of readscount in each bin at Input bam files}
}
\details{
This function is supposed to be called through the \code{Trumpet_report}
  wrapper. It is better not to call this function directly.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
