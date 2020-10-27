#include <Rcpp.h>
#include "mean_isfinite.h"

#ifndef _loki_kinmatrix_
#define _loki_kinmatrix_

template<typename scalar_t>
using VECTOR4 = Eigen::Matrix<scalar_t, 4, 1>;

template<typename scalar_t>
using VECTOR2 = Eigen::Matrix<scalar_t, 2, 1>;

template<typename scalar_t, class COEFF>
class KinMatrix {
  int size;
  public:
  std::vector<typename COEFF::template diagCoeff<scalar_t>> diagonale;
  std::vector<typename COEFF::template offDiagCoeff<scalar_t>> offDiagonale;
  KinMatrix(unsigned int n) : size(n), diagonale(n), offDiagonale( (n*(n-1))/2 ) {};
  void updateAdd(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2);
  void updateDom(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2);
  Rcpp::NumericMatrix getRawMatrix();
  Rcpp::NumericMatrix getInterceptMatrix();
};



// prend un SNP et update toute la matrice
template <typename scalar_t, class COEFF>
void KinMatrix<scalar_t, COEFF>::updateAdd(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2) {
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
  // scalar_t u1 = alpha*(1 - mu);
  // scalar_t u2 = alpha*(2 - mu);

  int k = 0;
  for(int i = 0; i < size; i++) {
    // scalar_t Xi = (1-P1[i]-P2[i])*u0 + P1[i]*u1 + P2[i]*u2:
    scalar_t Xi = u0 + P1[i]*alpha + P2[i]*2*alpha;
    scalar_t vi = P1[i]*(1 - P1[i]) + 4*P2[i]*(1 - P2[i]) - 4*P1[i]*P2[i];
    diagonale[i].update(Xi*Xi, vi);
    for(int j = i+1; j < size; j++) {
      scalar_t Xj = u0 + P1[j]*alpha + P2[j]*2*alpha;
      scalar_t vj = P1[j]*(1 - P1[j]) + 4*P2[j]*(1 - P2[j]) - 4*P1[j]*P2[j];
      offDiagonale[k++].update(Xi*Xj, vi, vj);
    }
  }
}

// idem pour le calcul de la matrice de dominance
template <typename scalar_t, class COEFF>
void KinMatrix<scalar_t, COEFF>::updateDom(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2) {
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
    scalar_t vi = P1[i]*(1 - P1[i]) + 4*P2[i]*(1 - P2[i]) - 4*P1[i]*P2[i];
    diagonale[i].update(Gamma*Xi*Xi, vi);
    for(int j = i+1; j < size; j++) {
      scalar_t Xj = (1-P1[j]-P2[j])*u0 + P1[j]*u1 + P2[j]*u2;
      scalar_t vj = P1[j]*(1 - P1[j]) + 4*P2[j]*(1 - P2[j]) - 4*P1[j]*P2[j];
      offDiagonale[k++].update(Gamma*Xi*Xj, vi, vj);
    }
  }
}


// renvoie la matrice symétrique des coeffs bruts...
template <typename scalar_t, class COEFF>
Rcpp::NumericMatrix KinMatrix<scalar_t, COEFF>::getRawMatrix() {
  Rcpp::NumericMatrix R = Rcpp::no_init_matrix(size, size);
  int k = 0;
  for(int i = 0; i < size; i++) {
    R(i,i) = diagonale[i].getPhi();
    for(int j = i+1; j < size; j++) {
      R(j,i) = offDiagonale[k++].getPhi();
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

// renvoie la matrice symétrique des intercepts
template <typename scalar_t, class COEFF>
Rcpp::NumericMatrix KinMatrix<scalar_t, COEFF>::getInterceptMatrix() {
  VECTOR2<scalar_t> C2;
  VECTOR4<scalar_t> C4;
  Rcpp::NumericMatrix R = Rcpp::no_init_matrix(size, size);
  int k = 0;
  /* diagonale[0].regress(C2);
  R(0,0) = C2(0);*/
  for(int i = 0; i < size; i++) {
    diagonale[i].regress(C2);
    R(i,i) = C2(0);
    for(int j = i+1; j < size; j++) {
      offDiagonale[k++].regress(C4);
      R(j,i) = C4(0);
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
