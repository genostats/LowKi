#include "milorGWAS/stringstream_lite.h"
#include "milorGWAS/token.h"
#include <cmath>
#ifndef VCFGP2PROBS
#define VCFGP2PROBS

// lit une chaîne du genre 0.998,0.002,0.000
// standardise en divisant par la somme
// et renvoie (p1, p2)
template<typename T>
std::pair<T, T> GP2probs(char * s) {
  if(*s == '.')
    return std::pair<T,T>(NAN,NAN);

  stringstream_lite ss(s, ',');
  T s0, s1, s2;
  if(ss.next_token())
    s0 = (T) atof(ss.token);
  else
    return std::pair<T,T>(NAN,NAN);

  if(ss.next_token())
    s1 = (T) atof(ss.token);
  else
    return std::pair<T,T>(NAN,NAN);

  if(ss.next_token())
    s2 = (T) atof(ss.token);
  else
    return std::pair<T,T>(NAN,NAN);

  s0 += s1+s2;

  return std::pair<T,T>(s1/s0, s2/s0);
}

#endif

