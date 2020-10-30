// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include "coeffRegression.h"
#include "coeffRegressionBis.h"
#include "coeff.h"
#include "KinMatrix.h"
using namespace Rcpp;

template<typename scalar_t, class C>
inline void fillKin(NumericMatrix P1, NumericMatrix P2, KinMatrix<scalar_t, C> & K, bool domi) {
  int n = P1.nrow();
  int m = P1.ncol();
  for(int j = 0; j < m; j++) {
    std::vector<scalar_t> p1(n);
    std::vector<scalar_t> p2(n);
    std::transform( &P1(0,j), &P1(0,j)+n, p1.begin(), [](double a) -> scalar_t { return (scalar_t) a; } );
    std::transform( &P2(0,j), &P2(0,j)+n, p2.begin(), [](double a) -> scalar_t { return (scalar_t) a; } );
//    std::copy( &P1(0,j), &P1(0,j)+n, p1.begin());
  //  std::copy( &P2(0,j), &P2(0,j)+n, p2.begin());
    if(domi)
      K.updateDom(p1, p2);
    else
      K.updateAdd(p1, p2);
  }
}

// (proof of principle)
// calcule la GRM et tous les coefficients de régression, en considant les colonnes de P1 et P2 une à une (SNP par SNP)
// (ce sont les matrices de probas des genotypes 1 et 2)
// [[Rcpp::export]]
NumericMatrix lowKincpp(NumericMatrix P1, NumericMatrix P2, bool adjust, bool domi, bool constr) {
  int n = P1.nrow();
  if(adjust) {
    if(constr) {
      KinMatrix<float, coeffRegressionBis> K(n);
      fillKin(P1, P2, K, domi);
      return K.getInterceptMatrix();
    } else {
      KinMatrix<float, coeffRegression> K(n);
      fillKin(P1, P2, K, domi);
      return K.getInterceptMatrix();
    }
  } else {
    KinMatrix<float, coeff> K(n);
    fillKin(P1, P2, K, domi);
    return K.getRawMatrix();
  }
}
