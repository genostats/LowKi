#include "milorGWAS/stringstream_lite.h"
#include "milorGWAS/token.h"
#include <cmath>
#ifndef PHREDSCORES
#define PHREDSCORES

// lit une chaîne du genre 0,3,35 
// 1 - fait la transformation (p0, p1, p2) = (10^0, 10^-0.3, 10^-3.5)
// 2 - standardise en divisant par la somme
// 3 - renvoie p1, p2
template<typename T>
std::pair<T, T> phred_P12(char * s) {
  if(*s == '.')
    return std::pair<T,T>(NAN,NAN);

  stringstream_lite ss(s, ',');
  int s0, s1, s2;
  if(ss.next_token())
    s0 = atoi(ss.token);
  else
    return std::pair<T,T>(NAN,NAN);

  if(ss.next_token())
    s1 = atoi(ss.token);
  else
    return std::pair<T,T>(NAN,NAN);

  if(ss.next_token())
    s2 = atoi(ss.token);
  else
    return std::pair<T,T>(NAN,NAN);

  T S, p1, p2;

  S = pow((T) 10.0, -0.1 * (T) s0);
  S += p1 = pow((T) 10.0, -0.1 * (T) s1);
  S += p2 = pow((T) 10.0, -0.1 * (T) s2);
  return std::pair<T,T>(p1/S, p2/S);
}

#endif

