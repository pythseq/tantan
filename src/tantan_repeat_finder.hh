// Copyright 2018 Martin C. Frith

// This class finds tandem repeats in a sequence, using a Viterbi
// algorithm.

// The input parameters are as in tantan.hh, with one difference: the
// substitutionMatrix should be a *log* likelihood ratio matrix:
// substitutionMatrix[x][y] = lambda * scoringMatrix[x][y].

// Usage: first call init.  Then call calcBestPathScore, which runs
// the Viterbi algorithm backwards from the end to the start of the
// sequence.  Finally, call nextState once per sequence letter, to get
// the "state" of each letter from start to end.

// state = 0: non-repeat.
// 0 < state <= maxRepeatOffset: tandem repeat with period = state.
// maxRepeatOffset < state < 2*maxRepeatOffset: insertion in repeat.

#ifndef TANTAN_REPEAT_FINDER_HH
#define TANTAN_REPEAT_FINDER_HH

#include <vector>

namespace tantan {

typedef unsigned char uchar;
typedef const double *const_double_ptr;

class RepeatFinder {
public:
  void init(int maxRepeatOffset,
	    const const_double_ptr *substitutionMatrix,
	    double repeatProb,
	    double repeatEndProb,
	    double repeatOffsetProbDecay,
	    double firstGapProb,
	    double otherGapProb);

  double calcBestPathScore(const uchar *seqBeg, const uchar *seqEnd);

  int nextState();

private:
  const const_double_ptr *substitutionMatrix;

  double b2b;
  double f2b;
  double g2g;
  double oneGapScore;
  double endGapScore;
  double f2f0;
  double f2f1;
  double f2f2;
  double b2fGrowth;
  double b2fLast;

  int maxRepeatOffset;
  int dpScoresPerLetter;

  std::vector<double> dpScores;
  double *scoresPtr;
  double *scoresEnd;
  double *checkpoint;

  const uchar *seqBeg;
  const uchar *seqEnd;
  const uchar *seqPtr;
  int state;

  void initializeBackwardScores();
  void calcBackwardTransitionScoresWithGaps();
  void calcEmissionScores();
  void calcScoresForOneSequencePosition();
  void makeCheckpoint();
  void redoCheckpoint();
  int offsetWithMaxScore() const;
  int deletionWithMaxScore() const;

  bool isNearSeqBeg() const {
    return seqPtr - seqBeg < maxRepeatOffset;
  }

  int maxOffsetInTheSequence() const {
    return isNearSeqBeg() ? (seqPtr - seqBeg) : maxRepeatOffset;
  }

  double scoreWithEmission(const double *matrixRow, int offset) const {
    return scoresPtr[offset] + matrixRow[seqPtr[-offset]];
  }
};

}

#endif
