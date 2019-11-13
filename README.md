[![Build Status](https://travis-ci.org/piotr-ole/kmer.svg?branch=master)](https://travis-ci.org/piotr-ole/kmer)
[![Codecov test coverage](https://codecov.io/gh/piotr-ole/kmer/branch/master/graph/badge.svg)](https://codecov.io/gh/piotr-ole/kmer?branch=master)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/piotr-ole/kmer?branch=master&svg=true)](https://ci.appveyor.com/project/piotr-ole/kmer)
# seqR <img src = "man/images/logo.png" align = "right" width="120"/>

Fast k-mer counting is crutial in working with peptides sequences. It may take part in buidling features for machine learining models for peptides prediction or as a part of sequence comparison software. SeqR provides fast k-mer counting implementation in R written with RCpp.

# Package installation and loading

Current verison is avaliable on Github, it can be build on Windows, Linux and Mac OS. For package installation from R Console the \code{devtools} package is needed. To install seqR and load it to R session the following instructions have to be executed.

```{r setup, eval=FALSE}
#install.packages('devtools') # if you dont have it installed
devtools::install_github("piotr-ole/seqR")
library(seqR)
```

# Examples

