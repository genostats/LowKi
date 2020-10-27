#ifndef _loki_coeff_offdiagcoeff_
#define _loki_coeff_offdiagcoeff_

template<typename scalar_t>
class coeff::offDiagCoeff {
  int n;
  scalar_t phi;
public:
  offDiagCoeff() : n(0), phi(0) {}

  scalar_t getPhi() {
    return phi / n;
  }
  
  void update(scalar_t phi_, scalar_t v_, scalar_t w_) {
    if(std::isfinite(phi_)) {
      n++;
      phi += phi_;
    }
  }
};

#endif

