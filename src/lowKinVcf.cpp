#include "vcf_reader.h"
#include "PL2probs.h"
#include "GP2probs.h"
#include "GenoProbas.h"
#include "coeffRegression.h"
#include "coeffRegressionBis.h"
#include "coeff.h"
#include "PartialKinMatrix.h"
#include "RawKinMatrix.h"
#include <Rcpp.h>

template<typename T, typename scalar_t, class C>
inline void fillKinVcf(vcf_reader<T> & VCF, PartialKinMatrix<scalar_t, C> & K, bool domi) {
  GenoProbas<scalar_t> probs( VCF.samples.size() );
  while(VCF.read_line(probs)) {
    if(domi) 
      K.updateDom(probs.P1, probs.P2);
    else
      K.updateAdd(probs.P1, probs.P2);
    probs.clear(); // remise à zero des vecteurs P1/P2 dans probs
  }
}

// !!!!!!!!! Attention le vecteur Index est utilisé pour des C++ index (premier élément = indice 0)
// [[Rcpp::export]]
NumericMatrix PartialKinVcf(std::string filename, IntegerVector Index, std::string field, bool adjust, bool domi, bool constr) {
  std::pair<float,float> (* CONVERT) (char *);
  if(field == "PL")
    CONVERT = PL2probs<float>;
  else if(field == "GP")
    CONVERT = GP2probs<float>;
  else
    stop("Unable to use field "+field);

  vcf_reader<std::pair<float,float>> VCF(filename, field, CONVERT);
  int n = VCF.samples.size();
  // check index...
  for(int i : Index) {
    if(i > n-1)
      stop("Index too large");
  }

  if(adjust) {
    if(constr) {
      PartialKinMatrix<float, coeffRegressionBis> K(Index);
      fillKinVcf(VCF, K, domi);
      return K.getInterceptMatrix();
    } else {
      PartialKinMatrix<float, coeffRegression> K(Index);
      fillKinVcf(VCF, K, domi);
      return K.getInterceptMatrix();
    }
  } else {
    PartialKinMatrix<float, coeff> K(Index);
    fillKinVcf(VCF, K, domi);
    return K.getRawMatrix();
  }
} 

// [[Rcpp::export]]
NumericMatrix RawKinVcf(std::string filename, std::string field, bool domi) {
  std::pair<float,float> (* CONVERT) (char *);
  if(field == "PL")
    CONVERT = PL2probs<float>;
  else if(field == "GP")
    CONVERT = GP2probs<float>;
  else
    stop("Unable to use field "+field);

  vcf_reader<std::pair<float,float>> VCF(filename, field, CONVERT);
  int n = VCF.samples.size();

  RawKinMatrix<float> K(n);
  GenoProbas<float> probs( VCF.samples.size() );

  while(VCF.read_line(probs)) {
    if(domi) 
      K.updateDom(probs.P1, probs.P2);
    else
      K.updateAdd(probs.P1, probs.P2);
    probs.clear(); // remise à zero des vecteurs P1/P2 dans probs
  }

  return K.getRawMatrix();
} 

