// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// aMTMsample
List aMTMsample(Function target, int N, int d, int K, arma::vec x0, arma::cube sig0, arma::mat mu0, arma::vec lam0, int adapt, int global, int scale, int local, int proposal, double accrate, double gamma, List parms, double beta);
RcppExport SEXP _aMTM_aMTMsample(SEXP targetSEXP, SEXP NSEXP, SEXP dSEXP, SEXP KSEXP, SEXP x0SEXP, SEXP sig0SEXP, SEXP mu0SEXP, SEXP lam0SEXP, SEXP adaptSEXP, SEXP globalSEXP, SEXP scaleSEXP, SEXP localSEXP, SEXP proposalSEXP, SEXP accrateSEXP, SEXP gammaSEXP, SEXP parmsSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type target(targetSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sig0(sig0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lam0(lam0SEXP);
    Rcpp::traits::input_parameter< int >::type adapt(adaptSEXP);
    Rcpp::traits::input_parameter< int >::type global(globalSEXP);
    Rcpp::traits::input_parameter< int >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< int >::type local(localSEXP);
    Rcpp::traits::input_parameter< int >::type proposal(proposalSEXP);
    Rcpp::traits::input_parameter< double >::type accrate(accrateSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< List >::type parms(parmsSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(aMTMsample(target, N, d, K, x0, sig0, mu0, lam0, adapt, global, scale, local, proposal, accrate, gamma, parms, beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_aMTM_aMTMsample", (DL_FUNC) &_aMTM_aMTMsample, 17},
    {NULL, NULL, 0}
};

RcppExport void R_init_aMTM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
