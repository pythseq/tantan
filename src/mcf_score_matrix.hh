// Copyright 2010 Martin C. Frith

// This class reads and writes scoring matrices in NCBI BLAST format.
// It can't be used directly: the idea is to copy it into a "fast"
// matrix, and then use that.

#ifndef MCF_SCORE_MATRIX_HH
#define MCF_SCORE_MATRIX_HH

#include <iosfwd>
#include <string>
#include <vector>

namespace mcf {

typedef unsigned char uchar;

class ScoreMatrix {
  // These preserve (upper/lower)case of letters labeling rows & columns.
  friend std::istream &operator>>(std::istream &s, ScoreMatrix &m);
  friend std::ostream &operator<<(std::ostream &s, const ScoreMatrix &m);

 public:
  static const char* blosum62;

  // Initializes the matrix with simple match/mismatch scores.
  void initMatchMismatch(int matchScore,  // usually > 0
                         int mismatchCost,  // usually > 0
                         const std::string &letters);  // case is preserved

  // Gets the minimum score in the matrix.
  int minScore() const;

  // Copies the ScoreMatrix into a "fast" matrix.
  // fastMatrixSize: the number of rows & columns in the fast matrix.
  // letterToIndex: mapping from letters, to indices in the fast
  // matrix.  If any index is >= fastMatrixSize, an error is thrown.
  // defaultScore: for cells in the fast matrix that lack any
  // corresponding cell in the ScoreMatrix.
  // isCaseSensitive: if false, then each cell in the ScoreMatrix will
  // be copied to uppercase and lowercase locations in the fast
  // matrix.
  // It is possible that multiple ScoreMatrix cells might map to one
  // fast matrix cell: there is currently no check for this, and later
  // ones will overwrite earlier ones.
  void makeFastMatrix(int **fastMatrix,
                      int fastMatrixSize,
                      const uchar *letterToIndex,
                      int defaultScore,
                      bool isCaseSensitive) const;

 private:
  std::vector<char> rowLetters;
  std::vector<char> colLetters;
  std::vector< std::vector<int> > scores;
};

}

#endif
