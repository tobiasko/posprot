% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.gcr.R
\name{get.MSstats_gcr}
\alias{get.MSstats_gcr}
\title{Extract ComparisonResult from a MSstats group comparision object}
\usage{
get.MSstats_gcr(x, label)
}
\arguments{
\item{x}{A group comparision result as generated by \code{\link[MSstats]{groupComparison}}}

\item{label}{Name of the comparision (character) to extract from x}
}
\value{
data.frame
}
\description{
This function extracts a named comparision result from an object returned by MSstats::groupComparison()
}
\details{
Acessor function for named $ComparisonResult. In case you would like to inpect which comparisions are available use \code{\link[posprot]{get.comparision.labels}}. An empty data.frame is returned if the requested comparision does not exist.
}
\examples{
n <- get.comparison.labels(sample.data)
n[1]
get.MSstats_gcr(sample.data, n[1])

}
