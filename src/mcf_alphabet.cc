// Copyright 2010 Martin C. Frith

#include "mcf_alphabet.hh"

#include <algorithm>  // fill_n
#include <cctype>
#include <stdexcept>

namespace mcf {

// consts have internal linkage, apparently.
// If I make these into static const class members, g++ 4.0.1 complains.
const uchar defaultNumber = 255;
const uchar defaultLetter = '!';

const char* Alphabet::dna = "ACGT";

// U=selenocysteine, O=pyrrolysine, *=stop?
const char* Alphabet::protein = "ACDEFGHIKLMNPQRSTVWY";

void Alphabet::fromString(const std::string &normalLetters) {
  std::fill_n(lettersToNumbers, capacity, defaultNumber);
  std::fill_n(numbersToLetters, capacity, defaultLetter);

  uchar number = 0;
  addLetters(normalLetters, number);
  size = number;

  addLetters("ABCDEFGHIJKLMNOPQRSTUVWXYZ", number);
  addLetters("abcdefghijklmnopqrstuvwxyz", number);
  addLetters("*", number);  // sometimes appears in protein sequences

  initCaseConversions();
}

void Alphabet::addLetters(const std::string &letters, uchar &number) {
  for (unsigned i = 0; i < letters.size(); ++i) {
    uchar letter = static_cast<uchar>(letters[i]);
    if (lettersToNumbers[letter] == defaultNumber) {
      lettersToNumbers[letter] = number;
      numbersToLetters[number] = letter;
      ++number;
    }
  }
}

void Alphabet::initCaseConversions() {
  for (unsigned i = 0; i < capacity; ++i) {
    uchar letter = numbersToLetters[i];
    numbersToUppercase[i] = lettersToNumbers[std::toupper(letter)];
    numbersToLowercase[i] = lettersToNumbers[std::tolower(letter)];
  }
}

void Alphabet::makeCaseInsensitive() {
  for (unsigned i = 0; i < capacity; ++i) {
    lettersToNumbers[i] = lettersToNumbers[std::toupper(i)];
  }
}

void Alphabet::encodeInPlace(uchar *beg, uchar *end) const {
  while (beg != end) {
    uchar letter = *beg;
    uchar number = lettersToNumbers[letter];
    if (number == defaultNumber) {
      throw std::runtime_error(std::string("bad symbol: ") + char(letter));
    }
    *beg = number;
    ++beg;
  }
}

void Alphabet::decodeInPlace(uchar *beg, uchar *end) const {
  while (beg != end) {
    uchar number = *beg;
    uchar letter = numbersToLetters[number];
    *beg = letter;
    ++beg;
  }
}

}
