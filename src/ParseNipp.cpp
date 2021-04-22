#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "ParseNipp.h"

NippEntry
parseNextGenus(std::ifstream& nippFile, const std::string & line)
{
  NippEntry entry;

  std::istringstream line_str(line);
  std::string desc;
  char next_char;

  line_str >> desc;
  #ifdef DEBUG
    assert( desc == "D=");
  #endif
  line_str >> entry.disc;

  line_str >> next_char;
  #ifdef DEBUG
    assert( next_char == ';');
  #endif

  line_str >> desc;
  #ifdef DEBUG
    assert( desc == "GENUS#");
  #endif
  line_str >> entry.genus;
  
  line_str >> next_char;
  #ifdef DEBUG
    assert(next_char == ';');
  #endif

  line_str >> desc;
  #ifdef DEBUG
    assert(desc == "MASS=");
  #endif
  line_str >> entry.mass[1];

  line_str >> next_char;
  #ifdef DEBUG
    assert(next_char == '/');
  #endif
  line_str >> entry.mass[2];

  line_str >> next_char;
  #ifdef DEBUG
    assert(next_char == ';');
  #endif

  line_str >> desc;
  #ifdef DEBUG
    assert(desc == "HASSE");
  #endif
  line_str >> desc;
  #ifdef DEBUG
    assert(desc == "SYMBOLS");
  #endif
  line_str >> desc;
  #ifdef DEBUG
    assert( desc == "ARE");
  #endif

  short int symb;
  while (line_str)
    {
      line_str >> symb;
      entry.HasseSymb.push_back(symb);
    }
  
  next_char = nippFile.peek();
  std::string next_line;
  // There might be trouble at the end of the file - !! TODO
  while (next_char != 'D')
    {
      LatticeRecord lattice;
      for (size_t i = 0; i < LatticeRecord::VecSize; i++)
	nippFile >> lattice.form[i];
      next_char = nippFile.get();
      #ifdef DEBUG
        assert(next_char == ';');
      #endif
      nippFile >> lattice.numAut;
      entry.lattices.push_back(lattice);
      next_char = nippFile.peek();
    }
  
  return entry;
}

std::vector<NippEntry>
parseDisc(const std::string & fname, const Z & disc)
{
  std::vector<NippEntry> genera;
  std::ifstream nippFile(fname);
  std::ostringstream disc_str, find_str;
  
  find_str << "D=";
  disc_str << disc;
  for (size_t i = 0; i < 5 - disc_str.str().size(); i++)
    find_str << " ";
  find_str << disc_str.str();

  if (nippFile.is_open())
    {
      std::string line;
      bool found = false;
      bool done = false;
      while (!done) && (std::getline(nippFile, line))
	{
	  size_t start = line.find(find_str, 0);
	  if (start != std::string::npos)
	    {
	     // found it
	     found = true;
	     NippEntry latGen = parseNextGenus(nippFile, line);
	     genera.push_back(latGen);
	    }
	  else
	    if (found) done = true; 
	}
    }
  else throw std::runtime_error("Unable to open nipp file.");
  nippFile.close();
  return genera;
}
