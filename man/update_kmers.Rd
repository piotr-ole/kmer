% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{update_kmers}
\alias{update_kmers}
\title{Update k-mers in a \code{hash map}}
\arguments{
\item{kmers}{a \code{hash map} reference representing the found k-mers (see details)}

\item{d}{\code{integer} vector representing gaps in k-mer}

\item{s}{\code{integer} vector representing an encoded sequence}

\item{kmer_hash}{\code{integer} representing computed hash of k-mer}

\item{kmer_begin_index}{\code{integer} representing the begin index of k-mer in \code{s}}

\item{num2str}{\code{hash map} representing encoding of sequence items between \code{integer} and \code{string}}

\item{kmer_decorator}{\code{function} that can add extra characters to \code{string} k-mer (for example position information)}
}
\description{
Update k-mers in a \code{hash map}
}
\details{
k-mers \code{hashmap} contains key-value pairs: key is an \code{integer} representing a hash of k-mer,
whereas the value represents a pair: (k-mer \code{string} for presentation, number of k-mer occurrences)
}
