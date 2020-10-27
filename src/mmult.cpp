#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix mmult(const NumericMatrix& m1) {
  NumericMatrix out(m1.nrow(),m1.nrow());
  NumericVector rm1, rm2;
  for (size_t i = 0; i < m1.nrow(); ++i) {
    rm1 = m1(i,_);
      for (size_t j = 0; j < m1.nrow(); ++j) {
        rm2 = m1(j,_);
        NumericVector w = rm1*rm2;
        NumericVector w2 = na_omit(w);
        out(i,j) = mean(w2);              
      }
    }
  return out;
}
