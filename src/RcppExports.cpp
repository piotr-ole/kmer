// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// get_window_length
int get_window_length(const Rcpp::IntegerVector& d);
RcppExport SEXP _kmer_get_window_length(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(get_window_length(d));
    return rcpp_result_gen;
END_RCPP
}
// get_hash
int get_hash(const std::vector<int>& s, const Rcpp::IntegerVector& d, int begin_index, bool pos);
RcppExport SEXP _kmer_get_hash(SEXP sSEXP, SEXP dSEXP, SEXP begin_indexSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type begin_index(begin_indexSEXP);
    Rcpp::traits::input_parameter< bool >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(get_hash(s, d, begin_index, pos));
    return rcpp_result_gen;
END_RCPP
}
// get_hash_for_word
int get_hash_for_word(const std::vector<int>& kmer);
RcppExport SEXP _kmer_get_hash_for_word(SEXP kmerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type kmer(kmerSEXP);
    rcpp_result_gen = Rcpp::wrap(get_hash_for_word(kmer));
    return rcpp_result_gen;
END_RCPP
}
// count_kmers_str
std::unordered_map<std::string, int> count_kmers_str(Rcpp::StringVector& s, Rcpp::IntegerVector& d, Rcpp::StringVector& alphabet, Rcpp::LogicalVector& pos);
RcppExport SEXP _kmer_count_kmers_str(SEXP sSEXP, SEXP dSEXP, SEXP alphabetSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type s(sSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type d(dSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type alphabet(alphabetSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector& >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(count_kmers_str(s, d, alphabet, pos));
    return rcpp_result_gen;
END_RCPP
}
// count_kmer_num
std::unordered_map<std::string, int> count_kmer_num(Rcpp::NumericVector& s, Rcpp::IntegerVector& d, Rcpp::NumericVector& alphabet, Rcpp::LogicalVector& pos);
RcppExport SEXP _kmer_count_kmer_num(SEXP sSEXP, SEXP dSEXP, SEXP alphabetSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type s(sSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type d(dSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type alphabet(alphabetSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector& >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(count_kmer_num(s, d, alphabet, pos));
    return rcpp_result_gen;
END_RCPP
}
// count_kmers_larger_than_one
std::unordered_map<std::string, int> count_kmers_larger_than_one(Rcpp::StringMatrix& m, Rcpp::IntegerVector& d, Rcpp::StringVector& alphabet, Rcpp::LogicalVector& pos);
RcppExport SEXP _kmer_count_kmers_larger_than_one(SEXP mSEXP, SEXP dSEXP, SEXP alphabetSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix& >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type d(dSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type alphabet(alphabetSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector& >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(count_kmers_larger_than_one(m, d, alphabet, pos));
    return rcpp_result_gen;
END_RCPP
}
// count_unigrams
std::unordered_map<std::string, int> count_unigrams(Rcpp::StringMatrix& m, Rcpp::StringVector& alphabet, Rcpp::LogicalVector& pos);
RcppExport SEXP _kmer_count_unigrams(SEXP mSEXP, SEXP alphabetSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix& >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type alphabet(alphabetSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector& >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(count_unigrams(m, alphabet, pos));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_kmer_get_window_length", (DL_FUNC) &_kmer_get_window_length, 1},
    {"_kmer_get_hash", (DL_FUNC) &_kmer_get_hash, 4},
    {"_kmer_get_hash_for_word", (DL_FUNC) &_kmer_get_hash_for_word, 1},
    {"_kmer_count_kmers_str", (DL_FUNC) &_kmer_count_kmers_str, 4},
    {"_kmer_count_kmer_num", (DL_FUNC) &_kmer_count_kmer_num, 4},
    {"_kmer_count_kmers_larger_than_one", (DL_FUNC) &_kmer_count_kmers_larger_than_one, 4},
    {"_kmer_count_unigrams", (DL_FUNC) &_kmer_count_unigrams, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_kmer(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
