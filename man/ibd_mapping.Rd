% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ibd_mapping.R
\name{ibd_mapping}
\alias{ibd_mapping}
\title{Conduct IBD mapping for a single local IBD matrix}
\usage{
ibd_mapping(
  glmmkin.randomvec.obj,
  Sigma = NULL,
  offset = NULL,
  ibd.ids,
  ibd.file
)
}
\arguments{
\item{glmmkin.randomvec.obj}{The glmmkin.randomvec object}

\item{Sigma}{An N by N sparse covariance matrix for the fitted glmmkin (default = NULL)}

\item{offset}{An N.randomvec by N.randomvec dense matrix computed as the cross-product of Sigma with glmmkin.randomvec matrix of random vectors on both sides (default = NULL)}

\item{ibd.ids}{The ID vector for row and column names of the local IBD matrix}

\item{ibd.file}{The input file names for the local IBD matrix}
}
\value{
A data.frame
\item{chr}{The chromosome}
\item{start}{The start position}
\item{end}{The end position}
\item{n.nonzero}{Number of non-zero entries in the upper triangle of the local IBD matrix}
\item{p.value}{Finite-sample p-value}
}
\description{
Conduct IBD mapping for a single local IBD matrix
}
