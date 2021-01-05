#include "phredScores.h"
#include <Rcpp.h>

//[[Rcpp::export]]
NumericVector bebopalula(std::string s) {
  std::pair<double, double> R = phred_P12<double>((char *) s.c_str());
  return NumericVector::create(1 - R.first - R.second, R.first, R.second);
}
