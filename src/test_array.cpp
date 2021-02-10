#include <Rcpp.h>
#include "multidimArray.h"
using namespace Rcpp;

// [[Rcpp::export]]
double essai_array(NumericVector A, IntegerVector I) {
  multidimArray<NumericVector, IntegerVector> x(A, A.attr("dim"));
  return x(I);
}
