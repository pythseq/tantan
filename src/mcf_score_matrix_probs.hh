// Copyright 2010 Martin C. Frith

// This class calculates the scale factor (lambda), and the letter
// probabilities, that are implicit in a scoring matrix.  The
// calculation might fail for various reasons, putting it into a
// "bad/undefined" state.

// If the score matrix is symmetric, then the two sets of letter
// probabilities should be identical.  With this code, they might
// differ minutely from exact identity.

#ifndef MCF_SCORE_MATRIX_PROBS_HH
#define MCF_SCORE_MATRIX_PROBS_HH

#include <vector>

namespace mcf {

typedef const int *const_int_ptr;

class ScoreMatrixProbs {
 public:
  ScoreMatrixProbs() { setBad(); }

  ScoreMatrixProbs(const const_int_ptr *scoreMatrix, unsigned alphabetSize)
  { init(scoreMatrix, alphabetSize); }

  void init(const const_int_ptr *scoreMatrix, unsigned alphabetSize);

  // Put us in the bad/undefined state.
  void setBad();

  // Are we in the bad/undefined state?
  bool isBad() const { return (lambda_ < 0); }

  // The scale factor.  In the bad/undefined state, it is negative.
  double lambda() const { return lambda_; }

  // The probabilities of letters corresponding to matrix rows (1st index).
  // In the bad/undefined state, it is an empty vector.
  const std::vector<double>& letterProbs1() const { return letterProbs1_; }

  // The probabilities of letters corresponding to matrix columns (2nd index).
  // In the bad/undefined state, it is an empty vector.
  const std::vector<double>& letterProbs2() const { return letterProbs2_; }

 private:
  double lambda_;
  std::vector<double> letterProbs1_;
  std::vector<double> letterProbs2_;
};

}

#endif
