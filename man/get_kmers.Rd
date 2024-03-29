% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{get_kmers}
\alias{get_kmers}
\title{Get k-mers}
\arguments{
\item{s}{input sequence}

\item{d}{\code{integer} vector representing gaps in k-mer}

\item{alphabet}{Rcpp sequence representing alphabet}

\item{val2str_converter}{\code{function} representing the conversion to string representation of an item}

\item{kmer_decorator}{\code{function} that can add extra characters in order to enhance the presentation of k-mer.}
}
\value{
\code{hash map} containing string representations of k-mers with their occurrence counts
}
\description{
Counts the occurrences of k-mers (the size of k-mer should be larger than one)
\code{B} - the (Rcpp) type of an input sequence
\code{S} - the (c++) type of an item of the sequence
}
