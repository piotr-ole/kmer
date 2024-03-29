% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{update_kmers_with_alphabet}
\alias{update_kmers_with_alphabet}
\title{Update k-mers with alphabet}
\arguments{
\item{kmers}{\code{hash map} reference that stores k-mers (see \link{update_kmers})}

\item{alphabet}{\code{integer} vector representing encoded alphabet}

\item{currentKmer}{\code{integer} vector representing the part of currently generated k-mer}

\item{k}{\code{integer} representing the number of k-mer items}

\item{num2str}{\code{hash map} representing the sequence encoding between \code{integer} and \code{string}}

\item{kmer_decorator}{\code{function} that can add extra characters for \code{string} representation of k-mer (for presentation reasons)}
}
\description{
Generates and add k-mers (based on the given alphabet) that do not exist in the \code{hash map}.
}
