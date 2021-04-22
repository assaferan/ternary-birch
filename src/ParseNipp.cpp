#include <fstream>
#include <string>
#include "ParseNipp.h"

std::vector<const QuadForm<Z,5> & > nippToForms(const NippEntry &)
{
  std::vector<const QuadForm<Z,5> & > ret;
  return ret;
}

NippEntry
parseNextGenus(std::ifstream& nippFile, const std::string & line)
{
  NippEntry entry;

  size_t start = 0;
  bool is_rec;
  std::istringstream line_str(line);
  std::string desc;
  char next_char;

  line_str >> desc;
  assert desc == "D=";
  line_str >> entry.disc;

  line_str >> next_char;
  assert next_char == ';';

  line_str >> desc;
  assert desc == "GENUS#";
  line_str >> entry.genus;
  
  line_str >> next_char;
  assert next_char == ';';

  line_str >> desc;
  assert desc == "MASS=";
  line_str >> entry.mass[1];

  line_str >> next_char;
  assert next_char == '/';
  line_str >> entry.mass[2];

  line_str >> next_char;
  assert next_char == ';';

  line_str >> desc;
  assert desc == "HASSE";
  line_str >> desc;
  assert desc == "SYMBOLS";
  line_str >> desc;
  assert desc == "ARE";

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
      for (size_t i = 0; i < lattice.form.size(); i++)
	nipp_file >> lattice.form[i];
      next_char = nippFile.get();
      assert next_char == ';';
      nipp_file >> lattice.numAut;
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

  if (nippFile.isopen())
    {
      size_t curLine = 0;
      std::string line;
      bool found = false;
      bool done = false;
      while (not done) and (std::getline(nippFile, line))
	{
	  curLine++;
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
