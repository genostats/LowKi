#include <Rcpp.h>
#include "mean_isfinite.h"

#ifndef _loki_rawkinmatrix_
#define _loki_rawkinmatrix_

/* same as KinMatrix but does not make the regression
 * mostly for allowing code similarity between functions
 * defined in lowKinVcf.cpp
 */

template<typename scalar_t>
class RawKinMatrix {
  int size;
  std::vector<scalar_t> Coeffs;
  std::vector<scalar_t> nbSnps; // pour chaque paire, le nb de snps pour lequel on a une estimation finie
  public:
  RawKinMatrix(unsigned int n) : size(n), Coeffs( (n*(n+1))/2, 0), nbSnps( (n*(n+1))/2, 0) {};
  void updateAdd(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2);
  void updateDom(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2);
  Rcpp::NumericMatrix getRawMatrix();
};


// prend un SNP et update toute la matrice
template <typename scalar_t>
void RawKinMatrix<scalar_t>::updateAdd(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2) {
  scalar_t mu1 = mean_isfinite(P1);
  scalar_t mu2 = mean_isfinite(P2);
  scalar_t mu0 = 1 - mu1 - mu2;
  // var genotypes in pop
  scalar_t Vad = mu1 + (4*mu0*mu2) - mu1*mu1;
  if(Vad == 0)
    return;
  // additive components
  scalar_t alpha = 1/std::sqrt(Vad);
  scalar_t mu = mu1 + 2*mu2;
  scalar_t u0 = alpha*(0 - mu);

  int k = 0;
  for(int i = 0; i < size; i++) {
     
    scalar_t Xi = u0 + P1[i]*alpha + P2[i]*2*alpha;
    if(std::isfinite(Xi)) {
      Coeffs[k] += Xi*Xi;
      nbSnps[k] += 1;
    }
    k++;
    if(std::isfinite(Xi)) {
      for(int j = i+1; j < size; j++) {
        scalar_t Xj = u0 + P1[j]*alpha + P2[j]*2*alpha;
        if(std::isfinite(Xj)) {
          Coeffs[k] += Xi*Xj;
          nbSnps[k] += 1;
        }
        k++;
      }
    } else {
      k += (size-i-1);
    }
  }
}

// idem pour le calcul de la matrice de dominance
template <typename scalar_t>
void RawKinMatrix<scalar_t>::updateDom(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2) {
  scalar_t mu1 = mean_isfinite(P1);
  scalar_t mu2 = mean_isfinite(P2);
  scalar_t mu0 = 1 - mu1 - mu2;

  scalar_t Gamma = 1/( mu0 + (4*mu0*mu2)/mu1 + mu2 );
  if(!std::isfinite(Gamma))
    return;
  
  // dominance components (up to sqrt(Gamma), reintroduced later)
  scalar_t u0 = sqrt(mu2/mu0);
  scalar_t u1 = -2*sqrt(mu0*mu2)/mu1;
  scalar_t u2 = 1/u0;
  
  int k = 0;
  for(int i = 0; i < size; i++) {
    scalar_t Xi = (1-P1[i]-P2[i])*u0 + P1[i]*u1 + P2[i]*u2;
    if(std::isfinite(Xi)) {
      Coeffs[k] += Gamma*Xi*Xi;
      nbSnps[k] += 1;
    }
    k++;
    if(std::isfinite(Xi)) {
      for(int j = i+1; j < size; j++) {
        scalar_t Xj = (1-P1[j]-P2[j])*u0 + P1[j]*u1 + P2[j]*u2;
        if(std::isfinite(Xj)) {
          Coeffs[k] += Gamma*Xi*Xj;
          nbSnps[k] += 1;
        }
        k++;
      }
    } else {
      k += (size-i-1);
    }
  }
}


// renvoie la matrice symétrique des coeffs (forcément Raw)
template <typename scalar_t>
Rcpp::NumericMatrix RawKinMatrix<scalar_t>::getRawMatrix() {
  Rcpp::NumericMatrix R = Rcpp::no_init_matrix(size, size);
  int k = 0;
  for(int i = 0; i < size; i++) {
    R(i,i) = (double) Coeffs[k] / (double) nbSnps[k];
    k++;
    for(int j = i+1; j < size; j++) {
      R(j,i) = (double) Coeffs[k] / (double) nbSnps[k];
      k++;
    }
  }
  // symétriser
  for(int i = 0; i < size; i++) {
    for(int j = 0; j < i; j++) {
      R(j,i) = R(i,j);
    }
  }
  return R;
}

#endif
