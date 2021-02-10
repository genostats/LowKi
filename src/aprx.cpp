#include <Rcpp.h>
#include "multidimArray.h"
#include "find_pos.h"
using namespace Rcpp;

// x et I ont longueur n 
// k = 0 au début, incrémente de 1 à chaque appel interne
// on considère les (n-k) derniers elts de x
// et les k premiers de I
// [[Rcpp::export]]
double aprx(NumericVector x, List L, NumericVector F_, IntegerVector ruleLeft, IntegerVector ruleRight, IntegerVector I, int k) {
  int n = x.size();
  if(k == n) {
    multidimArray<NumericVector, IntegerVector> F(F_, F_.attr("dim"));
    return F(I);
  }
  double xx = x[k];
  NumericVector X = as<NumericVector>(L[k]);
  int posxx = find_pos<NumericVector>(xx, X);
  int m = X.size();
  // est-ce qu'on est sorti de la grille à gauche ?
  if(posxx == -1) {
    switch(ruleLeft[k]) {
      case 1:
        return NA_REAL;
      case 2:
        I[k] = 0;
        return aprx(x, L, F_, ruleLeft, ruleRight, I, k+1);
      default:
        posxx = 0;
    }
  } else if(posxx == m-1) { // à drte
    switch(ruleRight[k]) { 
      case 1:
        return NA_REAL;
      case 2:
        I[k] = m-1;
        return aprx(x, L, F_, ruleLeft, ruleRight, I, k+1);
      default:
        posxx--;
     }
  }
  double x1 = X[posxx];
  double x2 = X[posxx+1];
  I[k] = posxx;
  double f1 = aprx(x, L, F_, ruleLeft, ruleRight, I, k+1);
  I[k] = posxx+1;
  double f2 = aprx(x, L, F_, ruleLeft, ruleRight, I, k+1);
  return ((x2 - xx)*f1 + (xx - x1)*f2)/(x2 - x1);
}
