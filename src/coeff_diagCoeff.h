#ifndef _loki_coeffregression_diag_coeff_
#define _loki_coeffregression diag_coeff_

template<typename scalar_t>
class coeff::diagCoeff {
  int n;
  scalar_t phi;
public:
  diagCoeff() : n(0), phi(0) { }

  scalar_t getPhi() {
    return phi / n;
  }
 
  void update(scalar_t phi_, scalar_t v_) {
    if(std::isfinite(phi_)) {
      n++;
      phi += phi_;
    }
  }
};

#endif 

