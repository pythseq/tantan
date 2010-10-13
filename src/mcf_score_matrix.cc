// Copyright 2010 Martin C. Frith

#include "mcf_score_matrix.hh"

#include <algorithm>  // min, fill_n
#include <cassert>
#include <cctype>
#include <iomanip>  // setw
#include <limits>
#include <sstream>
#include <stdexcept>

namespace mcf {

typedef std::runtime_error Error;

// copied from NCBI BLAST:
const char* ScoreMatrix::blosum62 = "\
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *\n\
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1 -1 -1 -4\n\
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1 -2  0 -1 -4\n\
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  4 -3  0 -1 -4\n\
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4 -3  1 -1 -4\n\
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -1 -3 -1 -4\n\
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0 -2  4 -1 -4\n\
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1 -3  4 -1 -4\n\
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -4 -2 -1 -4\n\
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0 -3  0 -1 -4\n\
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3  3 -3 -1 -4\n\
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4  3 -3 -1 -4\n\
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0 -3  1 -1 -4\n\
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3  2 -1 -1 -4\n\
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3  0 -3 -1 -4\n\
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -3 -1 -1 -4\n\
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0 -2  0 -1 -4\n\
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1 -1 -1 -4\n\
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -2 -2 -1 -4\n\
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -1 -2 -1 -4\n\
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3  2 -2 -1 -4\n\
B -2 -1  4  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4 -3  0 -1 -4\n\
J -1 -2 -3 -3 -1 -2 -3 -4 -3  3  3 -3  2  0 -3 -2 -1 -2 -1  2 -3  3 -3 -1 -4\n\
Z -1  0  0  1 -3  4  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -2 -2 -2  0 -3  4 -1 -4\n\
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\n\
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n\
";

void ScoreMatrix::initMatchMismatch(int matchScore,
                                    int mismatchCost,
                                    const std::string &letters) {
  rowLetters.assign(letters.begin(), letters.end());
  colLetters.assign(letters.begin(), letters.end());

  unsigned size = letters.size();
  assert(size == letters.size());

  scores.resize(size);

  for (unsigned i = 0; i < size; ++i) {
    scores[i].assign(size, -mismatchCost);
    scores[i][i] = matchScore;
  }
}

int ScoreMatrix::minScore() const {
  int m = std::numeric_limits<int>::max();
  for (unsigned i = 0; i < scores.size(); ++i)
    for (unsigned j = 0; j < scores[i].size(); ++j)
      m = std::min(m, scores[i][j]);
  return m;
}

static uchar indexOrDie(uchar letter, const uchar *letterToIndex, int tooBig) {
  uchar x = letterToIndex[letter];
  if (x >= tooBig)
    throw Error(std::string("bad letter in score matrix: ") +
                static_cast<char>(letter));
  return x;
}

void ScoreMatrix::makeFastMatrix(int **fastMatrix,
                                 int fastMatrixSize,
                                 const uchar *letterToIndex,
                                 int defaultScore,
                                 bool isCaseSensitive) const {
  for (int i = 0; i < fastMatrixSize; ++i)
    std::fill_n(fastMatrix[i], fastMatrixSize, defaultScore);

  for (unsigned i = 0; i < rowLetters.size(); ++i) {
    for (unsigned j = 0; j < colLetters.size(); ++j) {
      uchar x = static_cast<uchar>(rowLetters[i]);
      uchar y = static_cast<uchar>(colLetters[j]);

      if (isCaseSensitive) {
        uchar a = indexOrDie(x, letterToIndex, fastMatrixSize);
        uchar b = indexOrDie(y, letterToIndex, fastMatrixSize);
        fastMatrix[a][b] = scores[i][j];
      } else {
        uchar a = indexOrDie(std::toupper(x), letterToIndex, fastMatrixSize);
        uchar b = indexOrDie(std::tolower(x), letterToIndex, fastMatrixSize);
        uchar c = indexOrDie(std::toupper(y), letterToIndex, fastMatrixSize);
        uchar d = indexOrDie(std::tolower(y), letterToIndex, fastMatrixSize);
        fastMatrix[a][c] = scores[i][j];
        fastMatrix[a][d] = scores[i][j];
        fastMatrix[b][c] = scores[i][j];
        fastMatrix[b][d] = scores[i][j];
      }
    }
  }
}

std::istream &operator>>(std::istream &s, ScoreMatrix &m) {
  std::vector<char> tmpRowLetters;
  std::vector<char> tmpColLetters;
  std::vector< std::vector<int> > tmpScores;
  std::string line;

  while (std::getline(s, line)) {
    std::istringstream iss(line);
    char c;

    if (!(iss >> c)) continue;  // skip blank lines
    if (c == '#') continue;  // skip comment lines

    if (tmpColLetters.empty()) {
      do {
        tmpColLetters.push_back(c);
      } while (iss >> c);
    } else {
      tmpRowLetters.push_back(c);
      tmpScores.resize(tmpScores.size() + 1);
      int score;
      while (iss >> score) {
        tmpScores.back().push_back(score);
      }
      if (tmpScores.back().size() != tmpColLetters.size())
        s.setstate(std::ios::failbit);
    }
  }

  if (s.eof() && !s.bad() && !tmpScores.empty()) {
    s.clear();
    tmpRowLetters.swap(m.rowLetters);
    tmpColLetters.swap(m.colLetters);
    tmpScores.swap(m.scores);
  }

  return s;
}

std::ostream &operator<<(std::ostream &s, const ScoreMatrix &m) {
  s << ' ';
  for (unsigned j = 0; j < m.colLetters.size(); ++j)
    s << ' ' << std::setw(2) << m.colLetters[j];
  s << '\n';

  for (unsigned i = 0; i < m.rowLetters.size(); ++i) {
    s << m.rowLetters[i];
    for (unsigned j = 0; j < m.colLetters.size(); ++j) {
      s << ' ' << std::setw(2) << m.scores[i][j];
    }
    s << '\n';
  }

  return s;
}

}
