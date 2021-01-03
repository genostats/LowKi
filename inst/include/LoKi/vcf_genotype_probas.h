#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gaston/gzstream.h"
#include "milorGWAS/token.h"
#include "milorGWAS/read_gen_line.h"
#include "milorGWAS/read_vcf_header.h"
#include "milorGWAS/read_vcf_line.h"

#ifndef VCF_GENO_PROBA
#define VCF_GENO_PROBA

using namespace Rcpp;

class vcf_genotype_probas {
public:
  std::string filename;
  igzstream in;
  bool PHRED; // si vrai = champ PL, phred score ; si faux, champ GP, float
  std::string field; // should be "PL" or "GP"
  std::string line;
  bool good;
  std::vector<std::string> samples; 

  vcf_genotype_probas(std::string file, bool PHRED_);
  vcf_genotype_probas(std::string file, bool PHRED_, std::string field_);
  ~vcf_genotype_probas();


private:  
  void start();

public:
  // probas est un objet avec un membre .push_back(T1)
  // string_to_probas est un pointeur vers une fonction de 'T1 string_to_probas(std::string)'
  template<typename T0, typename T1> 
  bool read_line(T0 probas, T1 (*) (std::string) string_to_probas,
                 std::string & snp_id, int & snp_pos, std::string & chr, std::string & A1, std::string & A2);

  bool read_line(std::string & snp_id, int & snp_pos, std::string & chr, std::string & A1, std::string & A2);

};

#endif


