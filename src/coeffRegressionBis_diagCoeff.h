#include <RcppEigen.h>
#include <iostream>
#include "coeffRegression.h"

#ifndef _loki_coeffregressionbis_diag_coeff_
#define _loki_coeffregressionbis diag_coeff_

template<typename scalar_t>
using MATRIX2 = Eigen::Matrix<scalar_t, 2, 2>;

template<typename scalar_t>
using VECTOR2 = Eigen::Matrix<scalar_t, 2, 1>;

template<typename scalar_t>
class coeffRegressionBis::diagCoeff : public coeffRegression::diagCoeff<scalar_t> {
 // copie pure et simple de l'autre classe !
};

#endif 

