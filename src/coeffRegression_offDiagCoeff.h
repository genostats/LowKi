#include <RcppEigen.h>
#include <iostream>

#ifndef _loki_coeffregression_offdiagcoeff_
#define _loki_coeffregression_offdiagcoeff_

template<typename scalar_t>
using MATRIX4 = Eigen::Matrix<scalar_t, 4, 4>;

template<typename scalar_t>
using VECTOR4 = Eigen::Matrix<scalar_t, 4, 1>;

template<typename scalar_t>
class coeffRegression::offDiagCoeff {
  int n;
  scalar_t phi, v_phi, w_phi, vw_phi;
  scalar_t v, w, vw, v2, v2w, w2, vw2, v2w2;
public:
  offDiagCoeff() : n(0), phi(0), v_phi(0), w_phi(0), vw_phi(0), 
  v(0), w(0), vw(0), v2(0), v2w(0), w2(0), vw2(0), v2w2(0) { }

  scalar_t getPhi() {
    return phi / n;
  }
  
  void update(scalar_t phi_, scalar_t v_, scalar_t w_) {
    if(std::isfinite(phi_) && std::isfinite(v_) && std::isfinite(w_)) {
      scalar_t vw_ = v_*w_;
      scalar_t v2_ = v_*v_;
      scalar_t w2_ = w_*w_;
    
      n++;
      phi += phi_;
      v_phi += v_*phi_;
      w_phi += w_*phi_;
      vw_phi += vw_*phi_;
    
      v += v_;
      w += w_;
      vw += vw_;
      v2 += v2_;
      v2w += v2_*w_;
      w2 += w2_;
      vw2 += v_*w2_;
      v2w2 += v2_*w2_;
    }
  }
  
  void regress(VECTOR4<scalar_t> & coeffs) {
    VECTOR4<scalar_t> XtY;
    MATRIX4<scalar_t> XtX;
    XtY << phi, v_phi, w_phi, vw_phi;
    XtX << (scalar_t) n,  v,   w,   vw,
                      v,  v2,  vw,  v2w,
                      w,  vw,  w2,  vw2,
                      vw, v2w, vw2, v2w2;
    
    coeffs = XtX.colPivHouseholderQr().solve(XtY);
  }
};

#endif

