
#ifndef _loki_mean_isfinite_
#define _loki_mean_isfinite_
template<typename scalar_t>
scalar_t mean_isfinite(std::vector<scalar_t> x) {
  int n = 0;
  scalar_t S = 0;
  for(auto & a : x) {
    if(std::isfinite(a)) {
      S += a;
      n++;
    }
  }
  return S/n;
}
#endif

