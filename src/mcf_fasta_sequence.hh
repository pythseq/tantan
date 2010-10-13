// Copyright 2010 Martin C. Frith

#ifndef MCF_FASTA_SEQUENCE_HH
#define MCF_FASTA_SEQUENCE_HH

#include <iosfwd>
#include <string>
#include <vector>

namespace mcf {

typedef unsigned char uchar;

struct FastaSequence {
  std::string title;
  std::vector<uchar> sequence;
};

std::istream &operator>>(std::istream &s, FastaSequence &f);

std::ostream &operator<<(std::ostream &s, const FastaSequence &f);

}

#endif
