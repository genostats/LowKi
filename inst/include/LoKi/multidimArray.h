#include <Rcpp.h>
using namespace Rcpp;

#ifndef _multidimArray_
#define _multidimArray_

// A read-only array that allows to use R arrays
// (could be extended to r/w array if needed...)
template<typename dbvec, typename intvec>
class multidimArray {
public:
  dbvec x;
  intvec dim;
  // x is a vector data 
  // with x.size() == \prod_i dim(i)
  multidimArray(dbvec x_,intvec dim_) : x(x_), dim(dim_) {
    int l = 1;
    for(int a : dim) l *= a;
    if(l != x.size())
      stop("This is bad");
  }
  // pour aller plus vite on ne verifie pas que
  // pour tout i, I[i] < dim[i]
  template<typename intvec2>
  double operator()(intvec2 I) {
    int k = I[0];
    int l = 1;
    for(int i = 1; i < I.size(); i++) {
      l *= dim[i-1];
      k += l*I[i];
    }
    return x[k];
  }
};

#endif

/* EXAMPLE
// [[Rcpp::export]]
double essai_array(NumericVector A, IntegerVector I) {
 multidimArray<NumericVector, IntegerVector> x(A, A.attr("dim"));
 return x(I);
}
*/
