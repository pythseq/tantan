// Copyright 2010 Martin C. Frith

// An "Alphabet" consists of two things:
// 1) A mapping from ASCII letters to small integers.
// 2) A distinction between "normal" and "abnormal" letters.  The
// letters that map to the first "size" integers are normal.  For
// example, for DNA, ACGT would be normal.

#ifndef MCF_ALPHABET_HH
#define MCF_ALPHABET_HH

#include <string>

namespace mcf {

typedef unsigned char uchar;

struct Alphabet {
  static const char* dna;
  static const char* protein;

  static const unsigned capacity = 256;

  // Initializes the Alphabet with the given normal letters.
  void fromString(const std::string &normalLetters);  // case is preserved

  // Makes lowercase letters map to the same numbers that uppercase
  // letters do.
  void makeCaseInsensitive();

  // Maps letters to numbers in the given range.
  void encodeInPlace(uchar *beg, uchar *end) const;

  // Maps numbers back to letters in the given range.
  void decodeInPlace(uchar *beg, uchar *end) const;

  unsigned size;  // number of "normal" letters in the alphabet

  uchar lettersToNumbers[capacity];
  uchar numbersToLetters[capacity];

  uchar numbersToUppercase[capacity];
  uchar numbersToLowercase[capacity];

 private:
  void addLetters(const std::string &letters, uchar &number);
  void initCaseConversions();
};

}

#endif
