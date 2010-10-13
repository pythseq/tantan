// Copyright 2010 Martin C. Frith

#include "mcf_score_matrix_probs.hh"

namespace lambda {
extern "C" {
#include "CA_code/lambda_calculator.h"
}
}

namespace mcf {

void ScoreMatrixProbs::setBad() {
  lambda_ = -1;
  letterProbs1_.clear();
  letterProbs2_.clear();
}

void ScoreMatrixProbs::init(const const_int_ptr *scoreMatrix,
                            unsigned alphabetSize) {
  // We need to pass the parameters as 1-based pointers, hence the +1s
  // and -1s.

  std::vector<double> cells(alphabetSize * alphabetSize + 1);
  std::vector<const double*> pointers(alphabetSize + 1);

  for (unsigned i = 0; i < alphabetSize; ++i) {
    pointers[i+1] = &cells[i * alphabetSize];

    for (unsigned j = 0; j < alphabetSize; ++j) {
      cells[i * alphabetSize + j + 1] = scoreMatrix[i][j];
    }
  }

  letterProbs1_.resize(alphabetSize);
  letterProbs2_.resize(alphabetSize);

  lambda_ = lambda::calculate_lambda(&pointers[0], alphabetSize,
                                     &letterProbs1_[0] - 1,
                                     &letterProbs2_[0] - 1);

  if (lambda_ < 0) setBad();
}

}
