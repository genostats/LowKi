// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// lowKincpp
NumericMatrix lowKincpp(NumericMatrix P1, NumericMatrix P2, bool adjust, bool domi, bool constr);
RcppExport SEXP _LowKi_lowKincpp(SEXP P1SEXP, SEXP P2SEXP, SEXP adjustSEXP, SEXP domiSEXP, SEXP constrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P2(P2SEXP);
    Rcpp::traits::input_parameter< bool >::type adjust(adjustSEXP);
    Rcpp::traits::input_parameter< bool >::type domi(domiSEXP);
    Rcpp::traits::input_parameter< bool >::type constr(constrSEXP);
    rcpp_result_gen = Rcpp::wrap(lowKincpp(P1, P2, adjust, domi, constr));
    return rcpp_result_gen;
END_RCPP
}
// aprx
double aprx(NumericVector x, List L, NumericVector F_, IntegerVector ruleLeft, IntegerVector ruleRight, IntegerVector I, int k);
RcppExport SEXP _LowKi_aprx(SEXP xSEXP, SEXP LSEXP, SEXP F_SEXP, SEXP ruleLeftSEXP, SEXP ruleRightSEXP, SEXP ISEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type F_(F_SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ruleLeft(ruleLeftSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ruleRight(ruleRightSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type I(ISEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(aprx(x, L, F_, ruleLeft, ruleRight, I, k));
    return rcpp_result_gen;
END_RCPP
}
// KinVcf
NumericMatrix KinVcf(std::string filename, std::string field, bool adjust, bool domi, bool constr);
RcppExport SEXP _LowKi_KinVcf(SEXP filenameSEXP, SEXP fieldSEXP, SEXP adjustSEXP, SEXP domiSEXP, SEXP constrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type field(fieldSEXP);
    Rcpp::traits::input_parameter< bool >::type adjust(adjustSEXP);
    Rcpp::traits::input_parameter< bool >::type domi(domiSEXP);
    Rcpp::traits::input_parameter< bool >::type constr(constrSEXP);
    rcpp_result_gen = Rcpp::wrap(KinVcf(filename, field, adjust, domi, constr));
    return rcpp_result_gen;
END_RCPP
}
// PartialKinVcf
NumericMatrix PartialKinVcf(std::string filename, IntegerVector Index, std::string field, bool adjust, bool domi, bool constr);
RcppExport SEXP _LowKi_PartialKinVcf(SEXP filenameSEXP, SEXP IndexSEXP, SEXP fieldSEXP, SEXP adjustSEXP, SEXP domiSEXP, SEXP constrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Index(IndexSEXP);
    Rcpp::traits::input_parameter< std::string >::type field(fieldSEXP);
    Rcpp::traits::input_parameter< bool >::type adjust(adjustSEXP);
    Rcpp::traits::input_parameter< bool >::type domi(domiSEXP);
    Rcpp::traits::input_parameter< bool >::type constr(constrSEXP);
    rcpp_result_gen = Rcpp::wrap(PartialKinVcf(filename, Index, field, adjust, domi, constr));
    return rcpp_result_gen;
END_RCPP
}
// RawKinVcf
NumericMatrix RawKinVcf(std::string filename, std::string field, bool domi);
RcppExport SEXP _LowKi_RawKinVcf(SEXP filenameSEXP, SEXP fieldSEXP, SEXP domiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type field(fieldSEXP);
    Rcpp::traits::input_parameter< bool >::type domi(domiSEXP);
    rcpp_result_gen = Rcpp::wrap(RawKinVcf(filename, field, domi));
    return rcpp_result_gen;
END_RCPP
}
// mmult
NumericMatrix mmult(const NumericMatrix& m1);
RcppExport SEXP _LowKi_mmult(SEXP m1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type m1(m1SEXP);
    rcpp_result_gen = Rcpp::wrap(mmult(m1));
    return rcpp_result_gen;
END_RCPP
}
// readVcfProbas
List readVcfProbas(std::string filename, std::string field);
RcppExport SEXP _LowKi_readVcfProbas(SEXP filenameSEXP, SEXP fieldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type field(fieldSEXP);
    rcpp_result_gen = Rcpp::wrap(readVcfProbas(filename, field));
    return rcpp_result_gen;
END_RCPP
}
// essai_array
double essai_array(NumericVector A, IntegerVector I);
RcppExport SEXP _LowKi_essai_array(SEXP ASEXP, SEXP ISEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type I(ISEXP);
    rcpp_result_gen = Rcpp::wrap(essai_array(A, I));
    return rcpp_result_gen;
END_RCPP
}
// test_vcf_reader
List test_vcf_reader(std::string filename);
RcppExport SEXP _LowKi_test_vcf_reader(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(test_vcf_reader(filename));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LowKi_lowKincpp", (DL_FUNC) &_LowKi_lowKincpp, 5},
    {"_LowKi_aprx", (DL_FUNC) &_LowKi_aprx, 7},
    {"_LowKi_KinVcf", (DL_FUNC) &_LowKi_KinVcf, 5},
    {"_LowKi_PartialKinVcf", (DL_FUNC) &_LowKi_PartialKinVcf, 6},
    {"_LowKi_RawKinVcf", (DL_FUNC) &_LowKi_RawKinVcf, 3},
    {"_LowKi_mmult", (DL_FUNC) &_LowKi_mmult, 1},
    {"_LowKi_readVcfProbas", (DL_FUNC) &_LowKi_readVcfProbas, 2},
    {"_LowKi_essai_array", (DL_FUNC) &_LowKi_essai_array, 2},
    {"_LowKi_test_vcf_reader", (DL_FUNC) &_LowKi_test_vcf_reader, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_LowKi(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
