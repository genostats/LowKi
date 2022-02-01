#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include "vcf_reader.h"
#include "intPair.h"
#include "af_likelihood.h"

using namespace Rcpp;

// [[Rcpp::export]]
List vcfAlleleFreq(std::string filename) {
  vcf_reader<std::pair<int,int>> VCF( filename, "AD", intPair<int> );
  std::vector<std::string> SNP_ID, AL1, AL2;
  std::vector<int> CHR, POS;
  std::vector<double> P;
  std::vector< std::pair<int,int> > data;
  int nb_ind(-1); 
  while( VCF.read_line(data) ) {
    // check for right line size
    int len = data.size();
    if(nb_ind < 0) {
      nb_ind = len;
    } else if(nb_ind != len) {
      Rcerr << "While reading SNP #" << POS.size()+2 << " with id = " << VCF.snp_id << "\n";
      Rcerr << "Read " << len << " data fields, instead of " << nb_ind << " on previous line(s)\n";
      stop("File format error");
    }

    SNP_ID.push_back(VCF.snp_id);
    POS.push_back(VCF.snp_pos);
    CHR.push_back(VCF.chr);
    AL1.push_back(VCF.A1);
    AL2.push_back(VCF.A2);

    // estimate allele freq from data
    af_likelihood<double,int> LIK(data);
    // good starting point for Newton method
    std::pair<int,int> R = accumulate(data.begin(), data.end(), std::pair<int, int>(0,0),
              [](auto & a, auto & b){return std::pair<int, int>(a.first + b.first, a.second + b.second);});
    double p = ((double) R.first + 0.001) / ((double) R.first + (double) R.second + 0.001);

    bool b = LIK.newton_max(p, 0., 1., 1e-5, 10, false, true);
    if(!b) p = LIK.Brent_fmax(0., 1., 1e-5); // if Newton failed...

    P.push_back(1-p); // alternative allele freq !
    // clear data 
    data.clear();
  }
  List L;
  L["snp.id"] = wrap(SNP_ID);
  L["chr"] = wrap(CHR);
  L["pos"] = wrap(POS);
  L["A1"] = wrap(AL1);
  L["A2"] = wrap(AL2);
  L["p"] = wrap(P);
  return L;
}
