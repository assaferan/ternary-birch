#ifndef __NIPP_PARSE_H_
#define __NIPP_PARSE_H_

#include "QuadForm.h"

struct LatticeRecord {
  Z form[15];
  size_t numAut;
};

struct NippEntry {
  Z disc;
  size_t genus;
  Z mass[2]; // This is a rational number, stored as two integers
  // could use bitfields here, but that would complicate the code
  std::vector<short int> HasseSymb;
  std::vector<LatticeRecord> lattices;
};

class ParseNipp
{
 public:
  static std::vector<const NippEntry & >
    parseDisc(const std::string & fname, const Z & disc);
  
};

#endif // __NIPP_PARSE_H_
