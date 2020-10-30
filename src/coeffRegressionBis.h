#ifndef _loki_coeffregressionbis_
#define _loki_coeffregressionbis_

class coeffRegressionBis {
public:
  template <typename scalar_t>
  class offDiagCoeff;

  template <typename scalar_t>
  class diagCoeff;
};
#include "coeffRegressionBis_offDiagCoeff.h"
#include "coeffRegressionBis_diagCoeff.h"

#endif
