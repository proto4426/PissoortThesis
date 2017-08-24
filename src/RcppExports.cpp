// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// any_nonpos
bool any_nonpos(const Rcpp::NumericVector& x);
RcppExport SEXP _PissoortThesis_any_nonpos(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(any_nonpos(x));
    return rcpp_result_gen;
END_RCPP
}
// cpp_gev_loglik
double cpp_gev_loglik(const NumericVector x, const NumericVector data);
RcppExport SEXP _PissoortThesis_cpp_gev_loglik(SEXP xSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_gev_loglik(x, data));
    return rcpp_result_gen;
END_RCPP
}
// gev_logpost
double gev_logpost(const NumericVector& x, const NumericVector data);
RcppExport SEXP _PissoortThesis_gev_logpost(SEXP xSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(gev_logpost(x, data));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_statioCpp
List gibbs_statioCpp(NumericVector start, NumericVector data, int iter, NumericVector propsd, bool verbose);
RcppExport SEXP _PissoortThesis_gibbs_statioCpp(SEXP startSEXP, SEXP dataSEXP, SEXP iterSEXP, SEXP propsdSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propsd(propsdSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_statioCpp(start, data, iter, propsd, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PissoortThesis_any_nonpos", (DL_FUNC) &_PissoortThesis_any_nonpos, 1},
    {"_PissoortThesis_cpp_gev_loglik", (DL_FUNC) &_PissoortThesis_cpp_gev_loglik, 2},
    {"_PissoortThesis_gev_logpost", (DL_FUNC) &_PissoortThesis_gev_logpost, 2},
    {"_PissoortThesis_gibbs_statioCpp", (DL_FUNC) &_PissoortThesis_gibbs_statioCpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_PissoortThesis(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
