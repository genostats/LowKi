#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gaston/gzstream.h"
#include "gaston/flip_strand.h"
#include "milorGWAS/token.h"
#include "milorGWAS/stringstream_lite.h"
#include "milorGWAS/chr_convert.h"
#include "token2.h"

#ifndef GASTONread_vcf_line
#define GASTONread_vcf_line
using namespace Rcpp;

// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT

inline void parse_vcf_line(std::string line, std::string & snp_id, int & snp_pos, int & chr, 
                           std::string & A1, std::string & A2, double & qual, std::string & filter, std::string & info) {
  stringstream_lite li(line, 9); // 9 = tab separated
  std::string format, chr_;
  if(!(li >> chr_ >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format)) {
    stop("VCF file format error");
  }
  chr = chr_to_int(chr_);
}


template<typename T0, typename T1>
void parse_vcf_line_field(std::string line, std::string & snp_id, int & snp_pos, int & chr, 
                          std::string & A1, std::string & A2, double & qual, std::string & filter, std::string & info, 
                          std::string & field, T0 & data, T1 (*str2data) (char *)) {

  stringstream_lite li(line, 9); // 9 = tab separated
  std::string format, chr_;
  if(!(li >> chr_ >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format)) {
    stop("VCF file format error");
  }
  chr = chr_to_int(chr_);

  int pos = token_position(format, field);
  if(pos < 0) stop("VCF error (No field '" + field + "' found)");
  
  std::string D;
  while(li >> D) {
    T1 val = token_at_position<std::string>(D, pos, str2data);
    data.push_back(val);
  }
  return;
}

#endif
