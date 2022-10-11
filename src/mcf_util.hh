// Copyright 2010 Martin C. Frith

#ifndef MCF_UTIL_HH
#define MCF_UTIL_HH

#include "mcf_zstream.hh"

#include <cassert>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace mcf {

// open an input file, but if the name is "-", just return cin
std::istream &openIn(const std::string &fileName, izstream &z);

// open an output file, but if the name is "-", just return cout
std::ostream &openOut(const std::string &fileName, std::ofstream &ofs);

template <typename T>
std::string stringify(const T &x) {
  std::ostringstream oss;
  oss << x;
  assert(oss);
  return oss.str();
}

template <typename T>
void unstringify(T& x, const std::string &s) {
  std::istringstream iss(s);
  if (!(iss >> x) || !(iss >> std::ws).eof())
    throw std::runtime_error("can't interpret: " + s);
}

template <typename T>
void unfilify(T& x, const std::string &fileName) {
  izstream z;
  std::istream &input = openIn(fileName, z);
  input >> x;
  if (!input) throw std::runtime_error("can't read file: " + fileName);
  // check for junk at end of file?
}

}

#endif
