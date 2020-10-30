#include <RcppEigen.h>
#include <iostream>

#ifndef _loki_coeffregressionbis_offdiagcoeff_
#define _loki_coeffregressionbis_offdiagcoeff_

template<typename scalar_t>
using MATRIX3 = Eigen::Matrix<scalar_t, 3, 3>;

template<typename scalar_t>
using VECTOR3 = Eigen::Matrix<scalar_t, 3, 1>;

template<typename scalar_t>
class coeffRegressionBis::offDiagCoeff {
  int n;
  scalar_t phi, s_phi, p_phi;
  scalar_t s, p, s2, sp, p2;
public:
  offDiagCoeff() : n(0), phi(0), s_phi(0), p_phi(0),
  s(0), p(0), s2(0), sp(0), p2(0) {}

  scalar_t getPhi() {
    return phi / n;
  }
  
  void update(scalar_t phi_, scalar_t v_, scalar_t w_) {
    if(std::isfinite(phi_) && std::isfinite(v_) && std::isfinite(w_)) {
      scalar_t s_ = v_ + w_;
      scalar_t p_ = v_ * w_;

    
      n++;
      phi += phi_;
      s_phi += s_*phi_;
      p_phi += p_*phi_;
      
      s += s_;
      p += p_;
      s2 += s_ * s_;
      sp += s_ * p_;
      p2 += p_ * p_;
    } 
  }
  
  void regress(VECTOR3<scalar_t> & coeffs) {
    VECTOR3<scalar_t> XtY;
    MATRIX3<scalar_t> XtX;
    XtY << phi, s_phi, p_phi;
    XtX << (scalar_t) n, s,  p,
                      s, s2, sp,
                      p, sp, p2;
    
    coeffs = XtX.colPivHouseholderQr().solve(XtY);
  }

  scalar_t getIntercept() {
    VECTOR3<scalar_t> coeffs;
    regress(coeffs);
    return coeffs[0];
  }
  
};

#endif

