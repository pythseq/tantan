// Copyright 2010 Martin C. Frith

#include "mcf_fasta_sequence.hh"

#include <cctype>  // isspace
//#include <iostream>  // for debugging
#include <istream>
#include <iterator>  // istreambuf_iterator, ostreambuf_iterator
#include <ostream>

namespace mcf {

std::istream &operator>>(std::istream &s, FastaSequence &f) {
  std::string title;
  std::vector<unsigned char> sequence;

  {
    char c;
    while (s >> c)
      if (c == '>') break;
    getline(s, title);
  }
  if (!s) return s;

  std::istreambuf_iterator<char> inpos(s);
  std::istreambuf_iterator<char> endpos;
  while (inpos != endpos) {
    char c = *inpos;
    if (c == '>') break;  // we have hit the next FASTA sequence
    uchar u = static_cast<uchar>(c);
    if (!std::isspace(u)) sequence.push_back(u);
    ++inpos;
  }

  f.title.swap(title);
  f.sequence.swap(sequence);

  return s;
}

std::ostream &operator<<(std::ostream &s, const FastaSequence &f) {
  s << '>' << f.title << '\n';

  const int lettersPerLine = 50;  // ?

  std::ostreambuf_iterator<char> o(s);
  std::vector<uchar>::const_iterator b = f.sequence.begin();
  std::vector<uchar>::const_iterator e = f.sequence.end();
  std::vector<uchar>::const_iterator i = b;

  while (i != e) {
    o = static_cast<char>(*i);
    ++i;
    if (i - b == lettersPerLine || i == e) {
      o = '\n';
      b = i;
    }
  }

  return s;
}

}
