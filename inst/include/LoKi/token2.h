#include "milorGWAS/stringstream_lite.h"
/*
template<typename T>
T token_at_position(std::string & s, int pos, T (*str2data) (char *) ) {
  std::istringstream ss(s);
  std::string token;
  for(int i = 0; i < pos && std::getline(ss, token, ':'); i++) {}
  std::getline(ss, token, ':');
  T r = str2data((char *) token.c_str());
  return r;
}
*/

template<typename T>
T token_at_position(char * s, int pos, T (*str2data) (char *) ) {
  stringstream_lite ss(s, ':');
  std::string token;
  for(int i = 0; i <= pos && ss.next_token(); i++) {}
  T r = str2data((char *) ss.token);
  return r;
}

