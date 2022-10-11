// Copyright 2010 Martin C. Frith

#include "mcf_util.hh"

#include <iostream>

namespace mcf {

std::istream &openIn(const std::string &fileName, izstream &z) {
  if (fileName == "-") return std::cin;
  z.open(fileName.c_str());
  if (!z) throw std::runtime_error("can't open file: " + fileName);
  return z;
}

std::ostream &openOut(const std::string &fileName, std::ofstream &ofs) {
  if (fileName == "-") return std::cout;
  ofs.open(fileName.c_str());
  if (!ofs) throw std::runtime_error("can't open file: " + fileName);
  return ofs;
}

}
