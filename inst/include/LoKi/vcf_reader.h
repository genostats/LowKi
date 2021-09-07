#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gaston/gzstream.h"
#include "milorGWAS/token.h"
#include "milorGWAS/stringstream_lite.h"
#include "milorGWAS/chr_convert.h"
#include "token2.h"

#ifndef VCF_READER
#define VCF_READER

using namespace Rcpp;

template<typename T>
class vcf_reader {
private:
  std::string filename;
  igzstream in;
  std::string field;
  // str2data est un pointeur vers une fonction (char *) -> T
  // c'est la fonction qui prend le champ du vcf et la transforme en données de type T
  T (*str2data) (char *);
  bool keepAll;
  // keep un vecteur de booléens (std::vector<bool> ou LogicalVector...) 
  // pour ne lire qu'une partie des samples. C'est l'utilisateur qui 
  // veille à ce que 'keep' soit assez long.
  std::vector<bool> keep;
  bool good;
  std::string line;
public:
  // after calling the constructor, these vectors will be filled with
  // informations from the file prologue
  std::vector<std::string> info_ids, format_ids, samples; 
  // after each call to read_line, these elements will contain the corresponding fields of each line
  // chr_ is the string corresponding to the chromosome, chr is the result of chr_to_int(chr_)
  std::string chr_, snp_id, A1, A2, filter, info, format;
  int snp_pos, chr;
  double qual; 


  // file = vcf file
  // fi = field to read
  // str2data = fonction qui prend le 'field' lu (char *) et le tranforme en données de type T
  // NOTE on pourrait avoir un constructeur avec une valeur par défaut pour str2data, notamment
  // quand on ne veut utiliser que read_line() (la forme sans arguments)...
  vcf_reader(std::string file, std::string fi, T (*str2data_) (char *)) : 
       filename(file), in( (char *) &filename[0u] ), field(fi), str2data(str2data_), 
       keepAll(true), keep(std::vector<bool>(0)) {
    readHeader();
  }

  template<typename boolVec>
  vcf_reader(std::string file, std::string fi, T (*str2data_) (char *), boolVec & keep_) : 
       filename(file), in( (char *) &filename[0u] ), field(fi), str2data(str2data_), 
       keepAll(false) {
    // copie du vecteur booléen
    for(bool a : keep_) keep.push_back(a);
    readHeader();
  }

  void readHeader() {
    if(!in.good())
      stop("Can't open file");
    if(!std::getline(in, line))
      stop("File is empty");

    if(line.substr(0,1) != "#")
      stop("Not a VCF file");

    // go through description informations
    while(std::getline(in, line)) {
      if(line.substr(0,1) != "#") stop("Bad VCF format");

      if(line.substr(0,2) == "##") {
        std::string id;
        if(line.substr(0,11) == "##INFO=<ID=") { // read info fields
          std::istringstream li(line);
          std::getline(li, id, ',');
          info_ids.push_back(id.substr(11));
        } 
        if(line.substr(0,13) == "##FORMAT=<ID=") { // read format fields
          std::istringstream li(line);
          std::getline(li, id, ',');
          format_ids.push_back(id.substr(13));
        }
      } else { // read samples ids
        stringstream_lite li(line, 9); // 9 = tab sep.
        std::string G;
        for(int i = 0; i < 9; i++) { // skip col names
          if(!(li >> G))
         stop("VCF file format error");
        }
        while(li >> G) { // sample names
          samples.push_back( G );
        }
        break; // fin du prologue
      }
    }
    // compare length of sample / keep
    if(!keepAll && (keep.size() != samples.size())) {
      stop("length of 'keep' argument and # of samples mismatch");
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


  // data un objet (de type T0 quelconque) avec un membre .push_back(T)
  template<typename T0> 
  bool read_line(T0 & data) {
    if(!good) return false;

    stringstream_lite li(line, 9); // 9 = tab separated
    if(!(li >> chr_ >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format)) {
      stop("VCF file format error");
    }
    chr = chr_to_int(chr_);

    int pos = token_position(format, field);
    if(pos < 0) stop("VCF error (No field '" + field + "' found)");

    // deux boucles un peu différentes selon si keepAll vrai ou non
    if(keepAll) {
      while(li.next_token()) {
        T val = token_at_position(li.token, pos, str2data);
        data.push_back(val);
      } 
    } else {
      int i(0);
      while(li.next_token()) {
        if(keep[i++]) {
          T val = token_at_position(li.token, pos, str2data);
          data.push_back(val);
        }
      } 
    }

    if(std::getline(in, line))
      good = true;
    else
      good = false;
    return true;
  }


  // --------------------------------------------------------------------
  bool read_line() {
    if(!good) return false;

    stringstream_lite li(line, 9); // 9 = tab separated
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


