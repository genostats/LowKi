#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gaston/gzstream.h"
// #include "milorGWAS/read_gen_line.h"
#include "milorGWAS/read_vcf_header.h"
// #include "milorGWAS/read_vcf_line.h"
#include "parse_vcf_line_field.h"

#include "milorGWAS/token.h"
#include "milorGWAS/stringstream_lite.h"
#include "milorGWAS/chr_convert.h"

#ifndef VCF_READER
#define VCF_READER

using namespace Rcpp;

template<typename T>
class vcf_reader {
public:
  std::string filename;
  igzstream in;
  std::string field;
  T (*str2data) (char *);

  std::string line;
  bool good;
  std::vector<std::string> samples; 

  vcf_reader(std::string file, std::string fi, T (*str2data_) (char *)) : 
       filename(file), in( (char *) &filename[0u] ), field(fi), str2data(str2data_) {
    if(!in.good())
      stop("Can't open file");
    if(!std::getline(in, line))
      stop("File is empty");

    if(line.substr(0,1) != "#")
      stop("Not a VCF file");

    // skip description informations
    while(std::getline(in, line)) {
      if(line.substr(0,1) != "#") stop("Bad VCF format");
      if(line.substr(0,2) != "##") {
        read_vcf_samples(line, samples);
        break; // fin
      }
    }
    // read first data line
    if(std::getline(in, line))
      good = true;
    else
      good = false;
  }

  ~vcf_reader() {
    in.close();
  }


  // data un objet avec un membre .push_back(T1)
  // str2data est un pointeur vers une fonction de 'T1 string_to_probas(char *)'
  template<typename T0> 
  bool read_line(std::string & snp_id, int & snp_pos, int & chr, std::string & A1, std::string & A2, 
                 double & qual, std::string & filter, std::string & info, T0 & data) {
    if(!good) return false;

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
      T val = token_at_position(D, pos, str2data);
      data.push_back(val);
    }

    if(std::getline(in, line))
      good = true;
    else
      good = false;
    return true;
  }

  bool read_line(std::string & snp_id, int & snp_pos, int & chr, std::string & A1, std::string & A2, 
                 double & qual, std::string & filter, std::string & info) {
    if(!good) return false;

    stringstream_lite li(line, 9); // 9 = tab separated
    std::string format, chr_;
    if(!(li >> chr_ >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format)) {
      stop("VCF file format error");
    }
    chr = chr_to_int(chr_);

    if(std::getline(in, line))
      good = true;
    else
      good = false;
    return true;
  }
};

#endif


