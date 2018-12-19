// Copyright 2018 Martin C. Frith

#include "tantan_repeat_finder.hh"

#include <algorithm>
#include <assert.h>
#include <limits.h>
#include <math.h>
//#include <iostream>  // for debugging

namespace tantan {

static double max3(double x, double y, double z) {
  return std::max(std::max(x, y), z);
}

static double myLog(double x) {
  return x > 0 ? log(x) : -HUGE_VAL;
}

static double firstRepeatOffsetProb(double probMult, int maxRepeatOffset) {
  if (probMult < 1 || probMult > 1) {
    return (1 - probMult) / (1 - pow(probMult, 1.0 * maxRepeatOffset));
  }
  return 1.0 / maxRepeatOffset;
}

static int numOfDpScoresPerLetter(int maxRepeatOffset, double endGapScore) {
  if (endGapScore > -HUGE_VAL) {
    assert(maxRepeatOffset <= INT_MAX / 2);
    return maxRepeatOffset * 2;
  }
  assert(maxRepeatOffset < INT_MAX);
  return maxRepeatOffset + 1;
}

static unsigned minStoredPositions(const uchar *beg, const uchar *end) {
  // We will do a dynamic programming algorithm along the sequence.
  // To save memory, we will keep only some DP values, and recalculate
  // the others later.  We keep the values in an array x of size s.
  // We will store DP initial values in x[0].  We will put
  // the DP values after the first s-1 sequence positions in x[1],
  // the DP values after the next s-2 positions in x[2],
  // the DP values after the next s-3 positions in x[3], etc.
  // This function returns the minimum possible value of s.
  unsigned t = 0;
  while (t < end - beg) {
    beg += t;
    ++t;
  }
  return t + 1;
}

void RepeatFinder::init(int maxRepeatOffset,
			const const_double_ptr *substitutionMatrix,
			double repeatProb,
			double repeatEndProb,
			double repeatOffsetProbDecay,
			double firstGapProb,
			double otherGapProb) {
  assert(maxRepeatOffset > 0);
  this->maxRepeatOffset = maxRepeatOffset;
  this->substitutionMatrix = substitutionMatrix;

  b2b = myLog(1 - repeatProb);
  f2b = myLog(repeatEndProb);
  g2g = myLog(otherGapProb);
  oneGapScore = myLog(firstGapProb * (1 - otherGapProb));
  endGapScore = myLog(firstGapProb * (maxRepeatOffset > 1));
  f2f0 = myLog(1 - repeatEndProb);
  f2f1 = myLog(1 - repeatEndProb - firstGapProb);
  f2f2 = myLog(1 - repeatEndProb - firstGapProb * 2);

  double x = 1 / repeatOffsetProbDecay;
  b2fGrowth = myLog(x);
  b2fLast = myLog(repeatProb * firstRepeatOffsetProb(x, maxRepeatOffset));

  dpScoresPerLetter = numOfDpScoresPerLetter(maxRepeatOffset, endGapScore);
}

void RepeatFinder::initializeBackwardScores() {
  scoresPtr[0] = b2b;
  std::fill_n(scoresPtr + 1, maxRepeatOffset, f2b);
  if (endGapScore > -HUGE_VAL) {
    std::fill_n(scoresPtr + 1 + maxRepeatOffset, maxRepeatOffset-1, -HUGE_VAL);
  }
}

void RepeatFinder::calcBackwardTransitionScoresWithGaps() {
  double toBackground = f2b + scoresPtr[0];
  double *foregroundPtr = scoresPtr + 1;
  double f = *foregroundPtr;
  double toForeground = f;

  double *insertionPtr = scoresPtr + 1 + maxRepeatOffset;
  double i = *insertionPtr;
  *foregroundPtr = max3(toBackground, f2f1 + f, i);
  double d = endGapScore + f;
  ++foregroundPtr;
  toForeground += b2fGrowth;

  while (foregroundPtr < scoresPtr + maxRepeatOffset) {
    f = *foregroundPtr;
    toForeground = std::max(toForeground, f);
    i = *(insertionPtr + 1);
    *foregroundPtr = max3(toBackground, f2f2 + f, std::max(i, d));
    double oneGapScore_f = oneGapScore + f;
    *insertionPtr = std::max(oneGapScore_f, g2g + i);
    d = std::max(oneGapScore_f, g2g + d);
    ++foregroundPtr;
    ++insertionPtr;
    toForeground += b2fGrowth;
  }

  f = *foregroundPtr;
  toForeground = std::max(toForeground, f);
  *foregroundPtr = max3(toBackground, f2f1 + f, d);
  *insertionPtr = endGapScore + f;

  scoresPtr[0] = std::max(b2b + scoresPtr[0], b2fLast + toForeground);
}

void RepeatFinder::calcBackwardTransitionScores() {
  if (endGapScore > -HUGE_VAL) return calcBackwardTransitionScoresWithGaps();

  double toBackground = f2b + scoresPtr[0];
  double toForeground = -HUGE_VAL;
  double *foregroundPtr = scoresPtr + 1;
  double *foregroundEnd = foregroundPtr + maxRepeatOffset;

  while (foregroundPtr < foregroundEnd) {
    toForeground += b2fGrowth;
    double f = *foregroundPtr;
    toForeground = std::max(toForeground, f);
    *foregroundPtr = std::max(toBackground, f2f0 + f);
    ++foregroundPtr;
  }

  scoresPtr[0] = std::max(b2b + scoresPtr[0], b2fLast + toForeground);
}

void RepeatFinder::calcEmissionScores() {
  const double *matrixRow = substitutionMatrix[*seqPtr];
  double *oldScores = scoresPtr - dpScoresPerLetter;
  int maxOffset = maxOffsetInTheSequence();
  int i = 1;

  scoresPtr[0] = oldScores[0];

  for (; i <= maxOffset; ++i) {
    scoresPtr[i] = oldScores[i] + matrixRow[seqPtr[-i]];
  }

  for (; i <= maxRepeatOffset; ++i) {
    scoresPtr[i] = -HUGE_VAL;
  }

  std::copy(oldScores + i, scoresPtr, scoresPtr + i);
}

void RepeatFinder::makeCheckpoint() {
  checkpoint += dpScoresPerLetter;
  std::copy(scoresPtr - dpScoresPerLetter, scoresPtr, checkpoint);
  scoresPtr = checkpoint + dpScoresPerLetter;
  assert(scoresPtr < scoresEnd);
}

void RepeatFinder::redoCheckpoint() {
  seqPtr += (scoresEnd - scoresPtr) / dpScoresPerLetter;
  while (scoresPtr < scoresEnd) {
    --seqPtr;
    calcScoresForOneSequencePosition();
    scoresPtr += dpScoresPerLetter;
  }
  scoresPtr -= dpScoresPerLetter;
  checkpoint -= dpScoresPerLetter;
}

double RepeatFinder::calcBestPathScore(const uchar *seqBeg,
				       const uchar *seqEnd) {
  this->seqBeg = seqBeg;
  this->seqEnd = seqEnd;

  unsigned long numOfStoredPositions = minStoredPositions(seqBeg, seqEnd);
  unsigned long numOfScores = numOfStoredPositions * dpScoresPerLetter;
  assert(numOfStoredPositions > 0);
  assert(numOfScores > 0);
  dpScores.resize(numOfScores);
  scoresPtr = &dpScores[0];
  scoresEnd = scoresPtr + numOfScores;
  checkpoint = scoresPtr;
  seqPtr = seqEnd;

  initializeBackwardScores();

  while (seqPtr > seqBeg) {
    --seqPtr;
    scoresPtr += dpScoresPerLetter;
    if (scoresPtr == scoresEnd) makeCheckpoint();
    calcScoresForOneSequencePosition();
  }

  state = 0;
  return scoresPtr[0];
}

int RepeatFinder::offsetWithMaxScore() const {
  const double *matrixRow = substitutionMatrix[*seqPtr];
  int maxOffset = maxOffsetInTheSequence();
  int bestOffset = 0;
  double toForeground = -HUGE_VAL;

  for (int i = 1; i <= maxOffset; ++i) {
    toForeground += b2fGrowth;
    double f = scoreWithEmission(matrixRow, i);
    if (f > toForeground) {
      toForeground = f;
      bestOffset = i;
    }
  }

  return bestOffset;
}

int RepeatFinder::deletionWithMaxScore() const {
  const double *matrixRow = substitutionMatrix[*seqPtr];
  int bestOffset = 1;
  double f = scoreWithEmission(matrixRow, 1);
  double d = endGapScore + f;

  for (int i = 2; i < state; ++i) {
    d += g2g;
    f = scoreWithEmission(matrixRow, i);
    if (oneGapScore + f > d) {
      d = oneGapScore + f;
      bestOffset = i;
    }
  }

  return bestOffset;
}

int RepeatFinder::nextState() {
  double maxScore = scoresPtr[state];
  if (scoresPtr == checkpoint) redoCheckpoint();
  scoresPtr -= dpScoresPerLetter;

  if (state == 0) {
    if (b2b + scoresPtr[0] < maxScore) state = offsetWithMaxScore();
  } else if (state <= maxRepeatOffset) {
    if (f2b + scoresPtr[0] >= maxScore) {
      state = 0;
    } else if (endGapScore > -HUGE_VAL) {
      double f = scoreWithEmission(substitutionMatrix[*seqPtr], state);
      if (state == 1) {
	if (f2f1 + f < maxScore) state += maxRepeatOffset;
      } else if (state == maxRepeatOffset) {
	if (f2f1 + f < maxScore) state = deletionWithMaxScore();
      } else if (f2f2 + f < maxScore) {
	if (scoresPtr[state + maxRepeatOffset] >= maxScore) {
	  state += maxRepeatOffset;
	} else {
	  state = deletionWithMaxScore();
	}
      }
    }
  } else {
    ++state;
    if (state == dpScoresPerLetter || g2g + scoresPtr[state] < maxScore) {
      state -= maxRepeatOffset;
    }
  }

  ++seqPtr;
  return state;
}

}
