// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// effective_lengths
NumericVector effective_lengths(NumericVector txlengths, NumericVector rdlengths);
RcppExport SEXP _miso_effective_lengths(SEXP txlengthsSEXP, SEXP rdlengthsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type txlengths(txlengthsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rdlengths(rdlengthsSEXP);
    rcpp_result_gen = Rcpp::wrap(effective_lengths(txlengths, rdlengths));
    return rcpp_result_gen;
END_RCPP
}
// em_count
List em_count(NumericMatrix txreads, NumericVector txlengths, NumericVector weights, int ntx, int nr, unsigned int maxit, double reltol, double abstol);
RcppExport SEXP _miso_em_count(SEXP txreadsSEXP, SEXP txlengthsSEXP, SEXP weightsSEXP, SEXP ntxSEXP, SEXP nrSEXP, SEXP maxitSEXP, SEXP reltolSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type txreads(txreadsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type txlengths(txlengthsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type ntx(ntxSEXP);
    Rcpp::traits::input_parameter< int >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type reltol(reltolSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(em_count(txreads, txlengths, weights, ntx, nr, maxit, reltol, abstol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_miso_effective_lengths", (DL_FUNC) &_miso_effective_lengths, 2},
    {"_miso_em_count", (DL_FUNC) &_miso_em_count, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_miso(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
