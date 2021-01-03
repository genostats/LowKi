#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include "vcf_reader.h"

using namespace Rcpp;

template<typename scalar_t>
inline scalar_t geno_conv(char * s) {
  int le = strlen(s);
  scalar_t g = 0;
  if(le == 3) { // cas diploide
    if(s[0] == '1') g++;
    if(s[2] == '1') g++;
    if(s[0] == '.' || s[2] == '.') g = 3; // missing value : NA
  } else if(le == 1) { // cas haploide
    if(s[0] == '1') g++;
    if(s[0] == '.') g = 3; // missing value : NA
  } else {
    g = 3;
  }
  return g;
}

//[[Rcpp::export]]
List test_vcf_reader(std::string filename) {
  vcf_reader<int> in( filename, "GT", geno_conv<int> );
  std::string snp_id, A1, A2, filter, info;
  double qual;
  int chr, snp_pos;
  std::vector<std::string> SNP_ID, AL1, AL2;
  std::vector<int> CHR, POS;
  std::vector<double> data;
  int last_len(0), nb_ind(-1);
  while( in.read_line(snp_id, snp_pos, chr, A1, A2, qual, filter, info, data) ) {
    // check for right number of datas
    int len = data.size();
    if(nb_ind < 0) {
      nb_ind = len - last_len;
    } else if(nb_ind != len - last_len) {
      Rcerr << "While reading SNP #" << POS.size()+1 << " with id = " << snp_id << "\n";
      Rcerr << "Read " << len - last_len << " datas, instead of " << nb_ind << " on previous line(s)\n";
      stop("File format error");
    }
    last_len = len;

    SNP_ID.push_back(snp_id);
    POS.push_back(snp_pos);
    CHR.push_back(chr);
    AL1.push_back(A1);
    AL2.push_back(A2);
  }
  List L;
  L["snp.id"] = wrap(SNP_ID);
  L["chr"] = wrap(CHR);
  L["pos"] = wrap(POS);
  L["A1"] = wrap(AL1);
  L["A2"] = wrap(AL2);
  if(in.samples.size()> 0) L["samples"] = wrap(in.samples); // VCF ou PES
  NumericVector dat = wrap(data);
  dat.attr("dim") = Dimension( data.size()/POS.size(), POS.size() );
  L["data"] = dat;
  return L;
}

