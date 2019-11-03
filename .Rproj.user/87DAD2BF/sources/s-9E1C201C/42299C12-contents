#' @title Amino acid description dataframe
#'
#' @export
amino_acids <- function() {
  data.frame(Name = c("Arginine", "Histidine", "Lysine",
                       "Aspartic Acid", "Glutamic Acid",
                       "Serine", "Threonine", "Aspargine",
                       "Glutamine", "Cysteine", "Selenocysteine",
                       "Glycine", "Proline", "Alanine", "Valine",
                       "Isoleucine", "Leucine", "Methionine", "Phenylalanine",
                       "Tyrosine", "Tryptophan"),
             Short = c("Arg", "His", "Lys", "Asp", "Glu", "Ser",
                       "Thr", "Asn", "Gln", "Cys", "Sec", "Gly",
                       "Pro", "Ala", "Val", "Ile", "Leu", "Met",
                       "Phe", "Tyr", "Trp"),
             Letter = c("R", "H", "K", "D", "E", "S", "T", "N",
                        "Q", "C", "U", "G", "P", "A", "V", "I",
                        "L", "M", "F", "Y", "W"))
}

#' @title Amino acid letters
#'
#' @export
amino_letters <- function() {
  as.character(amino_acids()[, "Letter"])
}

#' @title Generating random alphabet of given size
#'
#' @param nsize  \code{integer} number of elements in the alphabet 
#'
#' @export
generate_random_alphabet <- function(nsize = 10) {
  set <- amino_letters()
  set[sample(length(set), size = 10)]
}

generate_random_seq <- function(nsize = 50) {
  set <- amino_letters()
  samp <- sample(set, nsize, replace = TRUE)
  paste(samp, collapse = '')
}


count_mers <- function(string) {
  l <- list()
  n <- nchar(string)
  for (j in 0:n) {
    for (i in 1:(n - j)) {
      l[j*n + i] <- substr(string, i, i+j)
    }
  }
  table(unlist(l))
}

count_k_mers <- function(string, k = 2) {
  l <- list()
  n <- nchar(string)
  for (i in 1:(n - k - 1)) {
      l[i] <- substr(string, i, i+k - 1)
    }
  table(unlist(l))
}

#'@title Count kmers
#'
#'@param string \code{character} with the sequence
#'@param k \code{integer} defines length of the subsequence to be counted
#'
#'@export
#'
#'@examples
#'count_k_mers_fast("AAPGAGAYY", 2)
count_k_mers_fast <- function(string, k) {
  reduce(map(string, k))
}


