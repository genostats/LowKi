#include "vcf_reader.h"
#include "PL2probs.h"
#include "GP2probs.h"
#include "GenoProbas.h"
#include <string>
#include <iostream>
#include <fstream>
#include <Rcpp.h>

//[[Rcpp::export]]
List readVcfProbas(std::string filename, std::string field) {
  std::pair<double,double> (* CONVERT) (char *);
  if(field == "PL")
    CONVERT = PL2probs<double>;
  else if(field == "GP")
    CONVERT = GP2probs<double>;
  else
    stop("Unable to use field "+field);

  vcf_reader<std::pair<double,double>> VCF(filename, field, CONVERT);
  std::vector<std::string> SNP_ID, AL1, AL2;
  std::vector<int> CHR, POS;
  std::vector<std::pair<double,double>> data;
  int last_len(0), nb_ind(-1);
  while( VCF.read_line(data) ) {
    // check for right number of datas
    int len = data.size();
    if(nb_ind < 0) {
      nb_ind = len - last_len;
    } else if(nb_ind != len - last_len) {
      Rcerr << "While reading SNP #" << POS.size()+1 << " with id = " << VCF.snp_id << "\n";
      Rcerr << "Read " << len - last_len << " datas, instead of " << nb_ind << " on previous line(s)\n";
      stop("File format error");
    }
    last_len = len;

    SNP_ID.push_back(VCF.snp_id);
    POS.push_back(VCF.snp_pos);
    CHR.push_back(VCF.chr);
    AL1.push_back(VCF.A1);
    AL2.push_back(VCF.A2);
  }

  // transformer le vecteur de std pair en deux matrices
  NumericMatrix P1(data.size()/POS.size(), POS.size());
  NumericMatrix P2(data.size()/POS.size(), POS.size());
  for(int i = 0; i < data.size(); i++) { 
    P1[i] = data[i].first;  // on accède aux éléments du vecteur sous jacent à P1...
    P2[i] = data[i].second;
  }

  List L;
  L["snp.id"] = wrap(SNP_ID);
  L["chr"] = wrap(CHR);
  L["pos"] = wrap(POS);
  L["A1"] = wrap(AL1);
  L["A2"] = wrap(AL2);
  L["samples"] = wrap(VCF.samples); 

  L["P1"] = P1;
  L["P2"] = P2;
  return L;
}

