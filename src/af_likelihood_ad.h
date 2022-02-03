#include "optimize_3.h"
#include <cmath>


#ifndef __loki_af_likelihood__
#define __loki_af_likelihood__

// pour l'optimisation de la vraisemblance basée sur AD (allele depth)
// T1 = float ou double, précision du calcul de la vraisemblance
// T2 = int ou char..., format des données
template<typename T1, typename T2>
class af_likelihood_ad : public fun_optim<T1> {
  std::vector<std::pair<T2, T2>> & data;
  T1 epsilon;
  public:
  af_likelihood_ad( std::vector<std::pair<T2, T2>> & data_, T1 eps = 0.001) : data(data_), epsilon(eps) {}

  // la log vraisemblance
  // l(p) = sum_{individus} log( p² Q0 + 2pq Q1 + q² Q2 ) avec Qi = P(data | genotype i)
  T1 f(T1 p) { 
    T1 lik = 0;
    for(std::pair<T2,T2> x : data) {
      // Qi = P( data | genotype i )
      T1 Q0 = std::pow(1-epsilon, (int) x.first) * std::pow(epsilon, (int) x.second); 
      T1 Q1 = std::pow(0.5, (int) x.first) * std::pow(0.5, (int) x.second);
      T1 Q2 = std::pow(epsilon, (int) x.first) * std::pow(1-epsilon, (int) x.second); 
      // p² * Q0 + 2p(1-p) Q1 + (1-p)² Q2 = p² (Q0 - 2Q1 + Q2) + 2p (Q1 - Q2) + Q2
      lik += log( Q2 + p*(2*(Q1 - Q2) + p*(Q0 - 2*Q1 + Q2)) );
    }
    return lik;
  }

  // ses dérivées
  void df_ddf(T1 p, T1 & df, T1 & ddf) {
    df = 0; 
    ddf = 0;
    for(std::pair<T2,T2> x : data) {
      // Qi = P( data | genotype i )
      T1 Q0 = std::pow(1-epsilon, (int) x.first) * std::pow(epsilon, (int) x.second); 
      T1 Q1 = std::pow(0.5, (int) x.first) * std::pow(0.5, (int) x.second);
      T1 Q2 = std::pow(epsilon, (int) x.first) * std::pow(1-epsilon, (int) x.second); 
      // p² * Q0 + 2p(1-p) Q1 + (1-p)² Q2 = p² (Q0 - 2Q1 + Q2) + 2p (Q1 - Q2) + Q2
      T1 R0 = Q0 - 2*Q1 + Q2;
      T1 R1 = Q1 - Q2; 
      T1 L = (Q2 + p*(2*R1 + p*R0));
      T1 dL = 2*(R1 + p*R0);
      T1 dlogL = dL/L;
      df += dlogL;
      ddf +=  2*R0/L - dlogL*dlogL;
    }
  }

};

#endif
