#include "vcf_reader.h"
#include "PL2probs.h"
#include "GenoProbas.h"
#include "coeffRegression.h"
#include "coeffRegressionBis.h"
#include "coeff.h"
#include "KinMatrix.h"
#include <Rcpp.h>

//[[Rcpp::export]]
NumericVector bebopalula(std::string s) {
  std::pair<double, double> R = PL2probs<double>((char *) s.c_str());
  return NumericVector::create(1 - R.first - R.second, R.first, R.second);
}


template<typename T, typename scalar_t, class C>
inline void fillKinVcf(vcf_reader<T> & VCF, KinMatrix<scalar_t, C> & K, bool domi) {
  GenoProbas<scalar_t> probs( VCF.samples.size() );
  while(VCF.read_line(probs)) {
    if(domi) 
      K.updateDom(probs.P1, probs.P2);
    else
      K.updateAdd(probs.P1, probs.P2);
    probs.clear(); // remise à zero des vecteurs P1/P2 dans probs
  }
}

// [[Rcpp::export]]
NumericMatrix lowKinVcf(std::string filename, std::string field, bool adjust, bool domi, bool constr) {
  std::pair<float,float> (* CONVERT) (char *);
  if(field == "PL")
    CONVERT = PL2probs<float>;
  else
    stop("Unable to use field "+field);

  vcf_reader<std::pair<float,float>> VCF( filename, "PL", CONVERT); //PL2probs<float> );
  int n = VCF.samples.size();
  if(adjust) {
    if(constr) {
      KinMatrix<float, coeffRegressionBis> K(n);
      fillKinVcf(VCF, K, domi);
      return K.getInterceptMatrix();
    } else {
      KinMatrix<float, coeffRegression> K(n);
      fillKinVcf(VCF, K, domi);
      return K.getInterceptMatrix();
    }
  } else {
    KinMatrix<float, coeff> K(n);
    fillKinVcf(VCF, K, domi);
    return K.getRawMatrix();
  }
} 

