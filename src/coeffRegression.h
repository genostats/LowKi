#ifndef _loki_coeffregression_
#define _loki_coeffregression_

class coeffRegression {
public:
  template <typename scalar_t>
  class offDiagCoeff;

  template <typename scalar_t>
  class diagCoeff;
};
#include "coeffRegression_offDiagCoeff.h"
#include "coeffRegression_diagCoeff.h"

#endif
