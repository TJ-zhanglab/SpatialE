% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getdelta.R
\name{getdelta}
\alias{getdelta}
\title{Calculating the delta values of all genes}
\usage{
getdelta(
  cluster_mean_exp_filtered,
  ExpMatrix,
  type = c("two sides", "high expression", "low expression")
)
}
\arguments{
\item{cluster_mean_exp_filtered}{a filtered clustered gene mean expression dataframe}

\item{ExpMatrix}{a matrix containing the expression values of all genes in all spots}

\item{type}{calculation type which to use. Available options are "two sides":calculating delta values in all clusters; "high expression":only calculating delta values of in clusters where mean expressions are above overall averages; "low expression":only calculating delta values of in clusters where mean expressions are below overall averages}
}
\value{
A dataframe containing delta values of all genes
}
\description{
Calculate the delta values (differences of one gene from the average expression in cluster to the overall average expression) of all genes, two-sides delta, high expression delta, low expression delta can be calculated
}
\examples{
delta <- getdelta(cluster_mean_exp_filtered, exp, type = "two sides")
}
