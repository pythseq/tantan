// Copyright 2010, 2011 Martin C. Frith

#include "mcf_fasta_sequence.hh"

#include <stddef.h>

#include <cctype>  // isspace
//#include <iostream>  // for debugging
#include <istream>
#include <iterator>  // istreambuf_iterator, ostreambuf_iterator
#include <ostream>

namespace mcf {

static void readSequence(std::istream &s, std::vector<uchar> &sequence,
                         char delimiter) {
  if (!s) return;
  std::istreambuf_iterator<char> inpos(s);
  std::istreambuf_iterator<char> endpos;
  while (inpos != endpos) {
    char c = *inpos;
    if (c == delimiter) break;
    uchar u = static_cast<uchar>(c);
    if (!std::isspace(u)) sequence.push_back(u);
    ++inpos;
  }
}

static void readQualityCodes(std::istream &s, std::vector<uchar> &qualityCodes,
                             std::vector<uchar>::size_type howMany) {
  if (!s) return;
  std::istreambuf_iterator<char> inpos(s);
  std::istreambuf_iterator<char> endpos;
  while (inpos != endpos) {
    if (qualityCodes.size() == howMany) break;
    char c = *inpos;
    uchar u = static_cast<uchar>(c);
    if (!std::isspace(u)) qualityCodes.push_back(u);
    ++inpos;
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
