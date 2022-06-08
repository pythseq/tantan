// Copyright 2010, 2011 Martin C. Frith

#include "mcf_fasta_sequence.hh"

#include <stddef.h>

//#include <iostream>  // for debugging
#include <istream>
#include <ostream>
#include <streambuf>

namespace mcf {

static void readSequence(std::istream &s, std::vector<uchar> &sequence,
                         char delimiter) {
  if (!s) return;
  std::streambuf *b = s.rdbuf();
  int c = b->sgetc();
  while (c != std::streambuf::traits_type::eof()) {
    if (c > ' ') {
      if (c == delimiter) break;
      sequence.push_back(c);
    }
    c = b->snextc();
  }
}

static void readQualityCodes(std::istream &s, std::vector<uchar> &qualityCodes,
                             std::vector<uchar>::size_type howMany) {
  if (!s) return;
  std::streambuf *b = s.rdbuf();
  while (howMany > 0) {
    int c = b->sbumpc();
    if (c == std::streambuf::traits_type::eof()) break;  // xxx ???
    if (c > ' ') {
      qualityCodes.push_back(c);
      --howMany;
    }
  }
}

std::istream &operator>>(std::istream &s, FastaSequence &f) {
  std::string title;
  std::vector<uchar> sequence;
  std::string secondTitle;
  std::vector<uchar> qualityCodes;

  char firstChar = '>';
  s >> firstChar;
  if (firstChar != '>' && firstChar != '@') s.setstate(std::ios::failbit);
  getline(s, title);

  if (firstChar == '>') {
    readSequence(s, sequence, '>');
  } else {
    readSequence(s, sequence, '+');
    char secondChar;
    s >> secondChar;
    getline(s, secondTitle);
    readQualityCodes(s, qualityCodes, sequence.size());
    // perhaps check whether we read enough quality codes
  }

  if (!s) return s;

  f.title.swap(title);
  f.sequence.swap(sequence);
  f.secondTitle.swap(secondTitle);
  f.qualityCodes.swap(qualityCodes);

  return s;
}

static void writeOneLine(std::ostream &s, const std::vector<uchar> &v) {
  std::streambuf *b = s.rdbuf();
  size_t size = v.size();
  if (size) b->sputn(reinterpret_cast<const char *>(&v[0]), size);
  b->sputc('\n');
}

static void writeMultiLines(std::ostream &s, const std::vector<uchar> &v) {
  size_t lettersPerLine = 50;  // ?
  std::streambuf *b = s.rdbuf();
  size_t size = v.size();
  for (size_t i = 0; i < size; i += lettersPerLine) {
    if (size - i < lettersPerLine) lettersPerLine = size - i;
    b->sputn(reinterpret_cast<const char *>(&v[i]), lettersPerLine);
    b->sputc('\n');
  }
}

std::ostream &operator<<(std::ostream &s, const FastaSequence &f) {
  if (f.qualityCodes.empty()) {
    s << '>' << f.title << '\n';
    writeMultiLines(s, f.sequence);
  } else {
    s << '@' << f.title << '\n';
    writeOneLine(s, f.sequence);
    s << '+' << f.secondTitle << '\n';
    writeOneLine(s, f.qualityCodes);
  }
  return s;
}

}
