#include <RcppEigen.h>
#include <iostream>

#ifndef _loki_coeffregression_diag_coeff_
#define _loki_coeffregression diag_coeff_

template<typename scalar_t>
using MATRIX2 = Eigen::Matrix<scalar_t, 2, 2>;

template<typename scalar_t>
using VECTOR2 = Eigen::Matrix<scalar_t, 2, 1>;

template<typename scalar_t>
class coeffRegression::diagCoeff {
  int n;
  scalar_t phi, v_phi, v, v2;
public:
  diagCoeff() : n(0), phi(0), v_phi(0), v(0), v2(0) { }

  scalar_t getPhi() {
    return phi / n;
  }

  // phi = estimateur à un SNP pour un individu avec lui meme
  // v = variance/fuziness en ce SNP
  // la régression est ~ 1 + v
 
  void update(scalar_t phi_, scalar_t v_) {
    if(std::isfinite(phi_) && std::isfinite(v_)) {
      n++;
      phi += phi_;
      v_phi += v_*phi_;
      v += v_;
      v2 += v_*v_;
    }
  }
  
  void regress(VECTOR2<scalar_t> & coeffs) {
    VECTOR2<scalar_t> XtY;
    MATRIX2<scalar_t> XtX;
    XtY << phi, v_phi;
    XtX << (scalar_t) n, v, 
                      v, v2;
//    Rcpp::Rcout << "n = " << n << "\n";
//    Rcpp::Rcout << XtY << "\n"; 
//    Rcpp::Rcout << XtX << "\n"; 
    coeffs = XtX.colPivHouseholderQr().solve(XtY);
//    Rcpp::Rcout << coeffs << "\n"; 
  }

  scalar_t getIntercept() {
    VECTOR2<scalar_t> coeffs;
    regress(coeffs);
    return coeffs[0];
  }
};

#endif 

