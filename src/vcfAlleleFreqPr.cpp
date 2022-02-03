#include <Rcpp.h>
#include "PL2probs.h"
#include "GP2probs.h"
#include "GenoProbas.h"
#include "mean_isfinite.h"
#include "vcf_reader.h"
#include <string>
#include <iostream>
#include <fstream>

using namespace Rcpp;

// [[Rcpp::export]]
List vcfAlleleFreqPr(std::string filename, std::string field) {
  std::pair<float,float> (* CONVERT) (char *);
  if(field == "PL")
    CONVERT = PL2probs<float>;
  else if(field == "GP")
    CONVERT = GP2probs<float>;
  else
    stop("Unable to use field "+field);

  vcf_reader<std::pair<float,float>> VCF(filename, field, CONVERT);

  std::vector<std::string> SNP_ID, AL1, AL2;
  std::vector<int> CHR, POS;
  std::vector<double> P;
  GenoProbas<float> data;
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
    float mu1 = mean_isfinite(data.P1);
    float mu2 = mean_isfinite(data.P2);
    P.push_back( 0.5*mu1 + mu2 );

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
