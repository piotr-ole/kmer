// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// countt_kmers
void countt_kmers(Rcpp::StringVector s, Rcpp::IntegerVector d, Rcpp::StringVector alphabet);
RcppExport SEXP _kmer_countt_kmers(SEXP sSEXP, SEXP dSEXP, SEXP alphabetSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type alphabet(alphabetSEXP);
    countt_kmers(s, d, alphabet);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_kmer_countt_kmers", (DL_FUNC) &_kmer_countt_kmers, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_kmer(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
