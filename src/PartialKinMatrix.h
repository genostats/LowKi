#include <Rcpp.h>
#include "mean_isfinite.h"

#ifndef _loki_partialkinmatrix_
#define _loki_partialkinmatrix_

/* computes the KinMatrix and does the regression of the coeffs on (vi vj) for a subset of individuals
 * the regression model is defined in the COEFF class.
 *
 * !! Still uses the whole probability vector to compute mu1 and mu2
 * An alternative would have been to give mu1 / mu2 to the update functions but
 * I don't think this would have been any better
 */

template<typename scalar_t, class COEFF>
class PartialKinMatrix {
  int size;
  std::vector<typename COEFF::template diagCoeff<scalar_t>> diagonale;
  std::vector<typename COEFF::template offDiagCoeff<scalar_t>> offDiagonale;
  std::vector<int> INDEX;
  public:
  template<typename intVec>
  PartialKinMatrix(intVec Index) : size(Index.size()), diagonale(size), offDiagonale( (size*(size-1))/2 ) {
    // copie du vecteur d'indices
    for(int a : Index) INDEX.push_back(a);
  }
  void updateAdd(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2);
  void updateDom(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2);
  void updateAdd(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2, scalar_t q);
  void updateDom(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2, scalar_t q);
  Rcpp::NumericMatrix getRawMatrix();
  Rcpp::NumericMatrix getInterceptMatrix();
};


/*******************************************
 *     les updates à la Loki               *
 *  Ils diffèrent des updates "classiques" *
 *  par leur prototype                     *
 *                                         *
 *******************************************/

// prend un SNP et update toute la matrice
template <typename scalar_t, class COEFF>
void PartialKinMatrix<scalar_t, COEFF>::updateAdd(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2) {
  /*** calcul des composantes additives à la Loki ***/
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
  for(int ii = 0; ii < size; ii++) {
    int i = INDEX[ii];
    //   P0 * u0 + P1 * u1 + P2 * u2 
    // = u0 + P1 * (u1 - u0) + P2 * (u2 - u0)
    // = u0 + P1 * alpha + P2 * 2 * alpha
    // = u0 + (P1 + 2 * P2) * alpha 
    scalar_t Xi = u0 + (P1[i] + 2*P2[i])*alpha;
    scalar_t vi = P1[i]*(1 - P1[i]) + 4*P2[i]*(1 - P2[i]) - 4*P1[i]*P2[i];
    diagonale[ii].update(Xi*Xi, vi); // attention on update diag[ii]
    for(int jj = ii+1; jj < size; jj++) {
      int j = INDEX[jj];
      scalar_t Xj = u0 + (P1[j] + 2*P2[j])*alpha;
      scalar_t vj = P1[j]*(1 - P1[j]) + 4*P2[j]*(1 - P2[j]) - 4*P1[j]*P2[j];
      offDiagonale[k++].update(Xi*Xj, vi, vj);
    }
  }
}

// idem pour le calcul de la matrice de dominance
template <typename scalar_t, class COEFF>
void PartialKinMatrix<scalar_t, COEFF>::updateDom(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2) {
  /*** calcul des composantes de dominance à la Loki ***/
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
  for(int ii = 0; ii < size; ii++) {
    int i = INDEX[ii];
    scalar_t Xi = (1-P1[i]-P2[i])*u0 + P1[i]*u1 + P2[i]*u2;
    scalar_t vi = P1[i]*(1 - P1[i]) + 4*P2[i]*(1 - P2[i]) - 4*P1[i]*P2[i];
    diagonale[ii].update(Gamma*Xi*Xi, vi); // attention on update diag[ii]
    for(int jj = ii+1; jj < size; jj++) {
      int j = INDEX[jj];
      scalar_t Xj = (1-P1[j]-P2[j])*u0 + P1[j]*u1 + P2[j]*u2;
      scalar_t vj = P1[j]*(1 - P1[j]) + 4*P2[j]*(1 - P2[j]) - 4*P1[j]*P2[j];
      offDiagonale[k++].update(Gamma*Xi*Xj, vi, vj);
    }
  }
}



/*******************************************
 *                                         *
 *     les updates "classiques"            *
 *  on leur passe aussi la fréquence du    *
 *  SNP                                    *
 *                                         *
 *******************************************/

// prend un SNP et update toute la matrice
template <typename scalar_t, class COEFF>
void PartialKinMatrix<scalar_t, COEFF>::updateAdd(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2, scalar_t q) {
  /*** calcul des composantes additives "classiques" ***/
  scalar_t Vad = 2*q*(1-q);
  if(Vad == 0)
    return;
  // additive components
  scalar_t alpha = 1/std::sqrt(Vad);
  scalar_t mu = 2*q;
  scalar_t u0 = alpha*(0 - mu);

  int k = 0;
  for(int ii = 0; ii < size; ii++) {
    int i = INDEX[ii];
    //   P0 * u0 + P1 * u1 + P2 * u2 
    // = u0 + P1 * (u1 - u0) + P2 * (u2 - u0)
    // = u0 + P1 * alpha + P2 * 2 * alpha
    // = u0 + (P1 + 2 * P2) * alpha 
    scalar_t Xi = u0 + (P1[i] + 2*P2[i])*alpha;
    scalar_t vi = P1[i]*(1 - P1[i]) + 4*P2[i]*(1 - P2[i]) - 4*P1[i]*P2[i];
    diagonale[ii].update(Xi*Xi, vi); // attention on update diag[ii]
    for(int jj = ii+1; jj < size; jj++) {
      int j = INDEX[jj];
      scalar_t Xj = u0 + (P1[j] + 2*P2[j])*alpha;
      scalar_t vj = P1[j]*(1 - P1[j]) + 4*P2[j]*(1 - P2[j]) - 4*P1[j]*P2[j];
      offDiagonale[k++].update(Xi*Xj, vi, vj);
    }
  }
}

// idem pour le calcul de la matrice de dominance
template <typename scalar_t, class COEFF>
void PartialKinMatrix<scalar_t, COEFF>::updateDom(const std::vector<scalar_t> & P1, const std::vector<scalar_t> & P2, scalar_t q) {
  /*** calcul des composantes de dominance "classique" ***/
  scalar_t u0 = q/(1-q);
  scalar_t u1 = -1.;
  scalar_t u2 = (1-q)/q;
  if(!std::isfinite(u0) || !std::isfinite(u1))
    return;
 
  int k = 0;
  for(int ii = 0; ii < size; ii++) {
    int i = INDEX[ii];
    scalar_t Xi = (1-P1[i]-P2[i])*u0 + P1[i]*u1 + P2[i]*u2;
    scalar_t vi = P1[i]*(1 - P1[i]) + 4*P2[i]*(1 - P2[i]) - 4*P1[i]*P2[i];
    diagonale[ii].update(Xi*Xi, vi); // attention on update diag[ii]
    for(int jj = ii+1; jj < size; jj++) {
      int j = INDEX[jj];
      scalar_t Xj = (1-P1[j]-P2[j])*u0 + P1[j]*u1 + P2[j]*u2;
      scalar_t vj = P1[j]*(1 - P1[j]) + 4*P2[j]*(1 - P2[j]) - 4*P1[j]*P2[j];
      offDiagonale[k++].update(Xi*Xj, vi, vj);
    }
  }
}





/*******************************************
 *                                         *
 *  getRawMatrix et getInterceptMatrix     *
 *                                         *
 *******************************************/


// renvoie la matrice symétrique des coeffs bruts...
template <typename scalar_t, class COEFF>
Rcpp::NumericMatrix PartialKinMatrix<scalar_t, COEFF>::getRawMatrix() {
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
Rcpp::NumericMatrix PartialKinMatrix<scalar_t, COEFF>::getInterceptMatrix() {
  Rcpp::NumericMatrix R = Rcpp::no_init_matrix(size, size);
  int k = 0;
  for(int i = 0; i < size; i++) {
    R(i,i) = diagonale[i].getIntercept();
    for(int j = i+1; j < size; j++) {
      R(j,i) = offDiagonale[k++].getIntercept();
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
