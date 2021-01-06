
// une classe pour stocker des triplets de probas génotypiques
// pas de membre x.P0, il faut prendre (1 - x.P1[i] - x.P2[i]) 
// Il y a des push_back (en particulier avec des pairs ou triplets
// comme argument, ce qui permet l'utilisation de la classe par vcf_reader)
// et les membres x.P12 et x.P012 pour récupérer des paires ou des
// triplets.

#ifndef GENOPROBAS
#define GENOPROBAS

template<typename T>
class GenoProbas {
public:
  std::vector<T> P1;
  std::vector<T> P2;

  GenoProbas() {};
  GenoProbas(int n) {
    P1.reserve(n);
    P2.reserve(n);
  }

  void push_back(T p1, T p2) {
    P1.push_back(p1);
    P2.push_back(p2);
  }
  void push_back(T p0, T p1, T p2) {
    p0 += p1+p2;
    P1.push_back(p1/p0);
    P2.push_back(p2/p0);
  }
  void push_back( std::pair<T, T> P12 ) {
    P1.push_back(P12.first);
    P2.push_back(P12.second);
  }
  void push_back( std::tuple<T, T, T> P012 ) {
    push_back( std::get<0>(P012), std::get<1>(P012), std::get<2>(P012) );
  }

  void clear() {
    P1.clear();
    P2.clear();
  }

  std::pair<T, T> P12(size_t i) {
    return std::make_pair(P1[i], P2[i]);
  }
  std::tuple<T, T, T> P012(size_t i) {
    T p1 = P1[i];
    T p2 = P2[i];
    return std::make_tuple(1.0-p1-p2, p1, p2);
  }

  size_t size() {
    return P1.size();
  }
};

#endif
