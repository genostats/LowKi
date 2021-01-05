#include "milorGWAS/stringstream_lite.h"
#ifndef TOKEN2
#define TOKEN2

template<typename T>
T token_at_position(char * s, int pos, T (*str2data) (char *) ) {
  stringstream_lite ss(s, ':');
  for(int i = 0; i <= pos && ss.next_token(); i++) {}
  T r = str2data((char *) ss.token);
  return r;
}

#endif

