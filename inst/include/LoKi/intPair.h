#include "milorGWAS/stringstream_lite.h"
#include "milorGWAS/token.h"
#ifndef _VCF_INTPAIR_
#define _VCF_INTPAIR_

// lit une chaîne du genre 0,3 -> std::pair<T,T>(0,3)
// les erreurs de lecture renvoient (0, 0) de façon silencieuse
template<typename T>
std::pair<T, T> intPair(char * s) {
  if(*s == '.')
    return std::pair<T,T>(0,0);

  stringstream_lite ss(s, ',');
  int s0, s1;
  if(ss.next_token())
    s0 = atoi(ss.token);
  else
    return std::pair<T,T>(0,0);

  if(ss.next_token())
    s1 = atoi(ss.token);
  else
    return std::pair<T,T>(0,0);

  return std::pair<T,T>((T) s0, (T) s1);
}

#endif

