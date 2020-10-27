#ifndef _loki_coeff_
#define _loki_coeff_

class coeff {
public:
  template <typename scalar_t>
  class offDiagCoeff;

  template <typename scalar_t>
  class diagCoeff;
};

#include "coeff_offDiagCoeff.h"
#include "coeff_diagCoeff.h"

#endif
