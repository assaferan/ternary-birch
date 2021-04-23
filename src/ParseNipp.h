#ifndef __PARSE_NIPP_H_
#define __PARSE_NIPP_H_

#include <string>
#include <vector>
#include "birch.h"

struct LatticeRecord {
  static const size_t n = 5;
  static const size_t VecSize = n*(n+1)/2;
  size_t form[VecSize];
  size_t numAut;
};

struct NippEntry {
  NippEntry() : HasseSymb(), lattices() {}
  size_t disc;
  size_t genus;
  size_t mass[2]; // This is a rational number, stored as two integers
  // could use bitfields here, but that would complicate the code
  std::vector<short int> HasseSymb;
  std::vector<LatticeRecord> lattices;
};

class ParseNipp
{
 public:
  static std::vector<NippEntry>
    parseDisc(const std::string & fname, const Z & disc);
  
};

#endif // __PARSE_NIPP_H_
