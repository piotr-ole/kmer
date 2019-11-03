library(testthat)

test_that('Correctness of alphabet', {
  expect_equal(kmer::amino_letters(), c("R", "H", "K", "D", "E", "S", "T", "N",
                                        "Q", "C", "U", "G", "P", "A", "V", "I",
                                        "L", "M", "F", "Y", "W"))
})