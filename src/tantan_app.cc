// Copyright 2010 Martin C. Frith

// Mask simple regions (low complexity & short-period tandem repeats)
// in biological sequences.

#include "mcf_alphabet.hh"
#include "mcf_fasta_sequence.hh"
#include "mcf_score_matrix.hh"
#include "mcf_tantan_options.hh"
#include "mcf_util.hh"
#include "tantan.hh"
#include "tantan_repeat_finder.hh"
#include "LambdaCalculator.hh"

#include <algorithm>  // copy, fill_n
#include <cassert>
#include <cmath>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <exception>  // exception
#include <fstream>
#include <iostream>
#include <new>  // bad_alloc
#include <string.h>

#define BEG(v) ((v).empty() ? 0 : &(v).front())
#define END(v) ((v).empty() ? 0 : &(v).back() + 1)

typedef std::runtime_error Error;

using namespace mcf;

// move this function to a reusable file?
bool isDubiousDna(const uchar *beg, const uchar *end) {
  int badLetters = 0;
  if (end - beg > 100) end = beg + 100;  // just check the first 100 letters
  while (beg < end) {
    if (!strchr("AaCcGgTtNnUu", *beg)) {
      ++badLetters;
      if (badLetters > 10) return true;
    }
    ++beg;
  }
  return false;
}

namespace {
TantanOptions options;
Alphabet alphabet;
tantan::RepeatFinder repeatFinder;

enum { scoreMatrixSize = 64 };
int fastMatrix[scoreMatrixSize][scoreMatrixSize];
int *fastMatrixPointers[scoreMatrixSize];
double probMatrix[scoreMatrixSize][scoreMatrixSize];
double *probMatrixPointers[scoreMatrixSize];

double firstGapProb;
double otherGapProb;

uchar hardMaskTable[Alphabet::capacity];
const uchar *maskTable;

std::vector<double> transitionCounts;
double transitionTotal;
}

void initAlphabet() {
  if (options.isProtein) alphabet.fromString(Alphabet::protein);
  else                   alphabet.fromString(Alphabet::dna);
}

void initScoresAndProbabilities() {
  std::copy(fastMatrix, fastMatrix + scoreMatrixSize, fastMatrixPointers);
  std::copy(probMatrix, probMatrix + scoreMatrixSize, probMatrixPointers);

  ScoreMatrix scoreMatrix;

  if (options.scoreMatrixFileName) {
    unfilify(scoreMatrix, options.scoreMatrixFileName);
  } else if (options.isProtein) {
    if (options.matchScore > 0 && options.mismatchCost > 0) {
      scoreMatrix.initMatchMismatch(options.matchScore, options.mismatchCost,
				    Alphabet::protein);
    } else {
      unstringify(scoreMatrix, ScoreMatrix::blosum62);
    }
  } else {
    scoreMatrix.initMatchMismatch(options.matchScore, options.mismatchCost,
				  "ACGTU");  // allow for RNA
  }

  scoreMatrix.makeFastMatrix(fastMatrixPointers, scoreMatrixSize,
                             alphabet.lettersToNumbers,
                             scoreMatrix.minScore(), false);

  cbrc::LambdaCalculator matCalc;
  matCalc.calculate(fastMatrixPointers, alphabet.size);
  if (matCalc.isBad())
    throw Error("can't calculate probabilities for this score matrix");
  double matrixLambda = matCalc.lambda();

  for (int i = 0; i < scoreMatrixSize; ++i) {
    for (int j = 0; j < scoreMatrixSize; ++j) {
      double x = matrixLambda * fastMatrix[i][j];
      if (options.outputType != options.repOut) x = std::exp(x);
      probMatrix[i][j] = x;
    }
  }

  if (options.gapExtensionCost > 0) {
    int firstGapCost = options.gapExistenceCost + options.gapExtensionCost;
    firstGapProb = std::exp(-matrixLambda * firstGapCost);
    otherGapProb = std::exp(-matrixLambda * options.gapExtensionCost);

    // the gap existence cost includes the gap ending cost:
    firstGapProb /= (1 - otherGapProb);

    // XXX check if firstGapProb is too high
  }

  repeatFinder.init(options.maxCycleLength, probMatrixPointers,
		    options.repeatProb, options.repeatEndProb,
		    options.repeatOffsetProbDecay, firstGapProb, otherGapProb);

  //std::cerr << "lambda: " << matrixLambda << "\n";
  //std::cerr << "firstGapProb: " << firstGapProb << "\n";
  //std::cerr << "otherGapProb: " << otherGapProb << "\n";
}

void initMaskTable() {
  uchar maskSymbol = static_cast<uchar>(options.maskSymbol);
  uchar maskNumber = alphabet.lettersToNumbers[maskSymbol];
  std::fill_n(hardMaskTable, Alphabet::capacity, maskNumber);

  if (maskSymbol == 0) maskTable = alphabet.numbersToLowercase;
  else                 maskTable = hardMaskTable;
}

std::string firstWord(const std::string &s) {
  std::string word;
  std::istringstream iss(s);
  iss >> word;
  return word;  // might be empty
}

void writeBedLine(const std::string &seqName, const float *origin,
                  const float *beg, const float *end, std::ostream &out) {
  out << seqName << '\t' << (beg - origin) << '\t' << (end - origin) << '\n';
}

void writeBed(const float *probBeg, const float *probEnd,
              const std::string &seqName, std::ostream &output) {
  if (seqName.empty()) throw Error("missing sequence name");
  const float *maskBeg = 0;  // pointer to start of masked tract
  for (const float *i = probBeg; i < probEnd; ++i) {
    if (*i >= options.minMaskProb) {  // this position is masked
      if (maskBeg == 0) maskBeg = i;
    } else {  // this position is not masked
      if (maskBeg) writeBedLine(seqName, probBeg, maskBeg, i, output);
      maskBeg = 0;
    }
  }
  if (maskBeg) writeBedLine(seqName, probBeg, maskBeg, probEnd, output);
}

void storeSequence(const uchar *beg, const uchar *end, std::string &out) {
  out.clear();
  for (const uchar *i = beg; i < end; ++i) {
    out.push_back(std::toupper(alphabet.numbersToLetters[*i]));
  }
}

struct RepeatUnit {
  const uchar *beg;
  int len;
};

bool less(const RepeatUnit &x, const RepeatUnit &y) {
  if (x.len != y.len) return x.len < y.len;
  int c = memcmp(x.beg, y.beg, x.len);
  if (c != 0) return c < 0;
  return x.beg > y.beg;
}

int mainLen(const std::vector<RepeatUnit> &repUnits) {
  int bestLen;
  size_t bestCount = 0;
  size_t count = 0;
  for (size_t i = 0; i < repUnits.size(); ++i) {
    if (i && repUnits[i].len > repUnits[i - 1].len) {
      count = 0;
    }
    ++count;
    if (count > bestCount) {
      bestCount = count;
      bestLen = repUnits[i].len;
    }
  }
  return bestLen;
}

const uchar *mainBeg(const std::vector<RepeatUnit> &repUnits, int len) {
  const uchar *bestBeg;
  size_t bestCount = 0;
  size_t count = 0;
  for (size_t i = 0; i < repUnits.size(); ++i) {
    if (repUnits[i].len != len) continue;
    if (count && memcmp(repUnits[i - 1].beg, repUnits[i].beg, len) != 0) {
      count = 0;
    }
    ++count;
    if (count < bestCount) continue;
    if (count > bestCount || repUnits[i].beg < bestBeg) {
      bestCount = count;
      bestBeg = repUnits[i].beg;
    }
  }
  return bestBeg;
}

void writeRepeat(const FastaSequence &f,
		 const uchar *repBeg, const uchar *repEnd,
		 const std::string &repText, std::vector<RepeatUnit> &repUnits,
		 const uchar *commaPos, int finalOffset) {
  double repeatCount = count(repText.begin(), repText.end(), ',');
  double copyNumber = repeatCount + (repEnd - commaPos) * 1.0 / finalOffset;
  if (copyNumber < options.minCopyNumber) return;

  sort(repUnits.begin(), repUnits.end(), less);
  int bestLen = mainLen(repUnits);
  const uchar *bestBeg = mainBeg(repUnits, bestLen);

  const uchar *beg = BEG(f.sequence);
  std::cout << firstWord(f.title) << '\t'
	    << (repBeg - beg) << '\t' << (repEnd - beg) << '\t'
	    << bestLen << '\t' << copyNumber << '\t';
  for (int i = 0; i < bestLen; ++i) {
    char c = std::toupper(alphabet.numbersToLetters[bestBeg[i]]);
    std::cout << c;
  }
  std::cout << '\t';
  std::cout << repText << '\n';
}

void findRepeatsInOneSequence(const FastaSequence &f) {
  const uchar *beg = BEG(f.sequence);
  const uchar *end = END(f.sequence);

  repeatFinder.calcBestPathScore(beg, end);

  std::vector<RepeatUnit> repUnits;
  std::string repText;
  const uchar *repBeg;
  const uchar *commaPos;
  int state = 0;

  for (const uchar *seqPtr = beg; seqPtr < end; ++seqPtr) {
    int newState = repeatFinder.nextState();

    if (newState == 0) {
      if (state > 0) {
	writeRepeat(f, repBeg, seqPtr, repText, repUnits, commaPos, state);
      }
    } else if (newState <= options.maxCycleLength) {
      if (state == 0) {
	repUnits.clear();
	repBeg = seqPtr - newState;
	storeSequence(repBeg, seqPtr, repText);
	commaPos = repBeg;
      } else if (state <= options.maxCycleLength) {
	for (int i = state; i > newState; --i) {
	  if (seqPtr - commaPos >= i) {
	    repText.push_back(',');
	    commaPos = seqPtr;
	  }
	  repText.push_back('-');
	}
      }
      RepeatUnit unit = {seqPtr - newState, newState};
      repUnits.push_back(unit);
      if (seqPtr - commaPos >= newState) {
	repText.push_back(',');
	commaPos = seqPtr;
      }
      repText.push_back(std::toupper(alphabet.numbersToLetters[*seqPtr]));
    } else {
      repText.push_back(std::tolower(alphabet.numbersToLetters[*seqPtr]));
    }

    state = newState;
  }

  if (state > 0) {
    writeRepeat(f, repBeg, end, repText, repUnits, commaPos, state);
  }
}

void processOneSequence(FastaSequence &f, std::ostream &output) {
  uchar *beg = BEG(f.sequence);
  uchar *end = END(f.sequence);

  alphabet.encodeInPlace(beg, end);

  if (options.outputType == options.maskOut) {
    tantan::maskSequences(beg, end, options.maxCycleLength,
                          probMatrixPointers,
                          options.repeatProb, options.repeatEndProb,
                          options.repeatOffsetProbDecay,
                          firstGapProb, otherGapProb,
                          options.minMaskProb, maskTable);
    alphabet.decodeInPlace(beg, end);
    output << f;
  } else if (options.outputType == options.countOut) {
    tantan::countTransitions(beg, end, options.maxCycleLength,
                             probMatrixPointers,
                             options.repeatProb, options.repeatEndProb,
                             options.repeatOffsetProbDecay,
                             firstGapProb, otherGapProb,
                             BEG(transitionCounts));
    double sequenceLength = static_cast<double>(f.sequence.size());
    transitionTotal += sequenceLength + 1;
  } else if (options.outputType == options.repOut) {
    findRepeatsInOneSequence(f);
  } else {
    std::vector<float> probabilities(end - beg);
    float *probBeg = BEG(probabilities);
    float *probEnd = END(probabilities);
    tantan::getProbabilities(beg, end, options.maxCycleLength,
                             probMatrixPointers,
                             options.repeatProb, options.repeatEndProb,
                             options.repeatOffsetProbDecay,
                             firstGapProb, otherGapProb, probBeg);
    if (options.outputType == options.probOut) {
      output << '>' << f.title << '\n';
      for (float *i = probBeg; i < probEnd; ++i)
        output << *i << '\n';
    } else {
      writeBed(probBeg, probEnd, firstWord(f.title), output);
    }
  }
}

void processOneFile(std::istream &input, std::ostream &output) {
  bool isFirstSequence = true;
  FastaSequence f;
  while (input >> f) {
    if (isFirstSequence && !options.isProtein &&
        isDubiousDna(BEG(f.sequence), END(f.sequence)))
      std::cerr << "tantan: that's some funny-lookin DNA\n";
    processOneSequence(f, output);
    isFirstSequence = false;
  }
}

void writeCounts(std::ostream &output) {
  double bg2bg = transitionCounts[0];

  output << "#period" << '\t' << "estimated number of tracts" << '\n';
  double repeatCountSum = 0;
  double weightedSum = 0;
  for (int i = 1; i <= options.maxCycleLength; ++i) {
    assert(transitionCounts[i] >= 0);
    output << i << '\t' << transitionCounts[i] << '\n';
    repeatCountSum += transitionCounts[i];
    weightedSum += i * transitionCounts[i];
  }

  output << "# estimated total number of repetitive tracts: "
         << repeatCountSum << '\n';

  output << "# estimated number of background-to-background transitions: "
         << bg2bg << '\n';

  output << "# total number of transitions (#letters + #sequences): "
         << transitionTotal << '\n';

  double probDecay = 1 - repeatCountSum / weightedSum;
  output << "# best-fit probability decay per period: " << probDecay << '\n';

  double repeatProb = repeatCountSum / (repeatCountSum + bg2bg);
  output << "# best-fit probability of a repeat starting per position: "
         << repeatProb << '\n';

  if (firstGapProb > 0) return;
  // The remaining calculations are correct only if there are no gaps.

  // If repeatProb and repeatEndProb are fixed to be the same value,
  // and there are no gaps:
  // repeatProb = 2 * repeatCountSum / transitionTotal;

  double fg2fg = transitionTotal - bg2bg - 2 * repeatCountSum;
  double repeatEndProb = repeatCountSum / (repeatCountSum + fg2fg);
  output << "# best-fit probability of a repeat ending per position: "
         << repeatEndProb << '\n';
}

int main(int argc, char **argv)
try {
  options.fromArgs(argc, argv);

  initAlphabet();
  initScoresAndProbabilities();
  initMaskTable();

  // do this after initMaskTable, so that the mask symbol can be lowercase:
  if (!options.isPreserveLowercase) alphabet.makeCaseInsensitive();

  if (options.outputType == options.countOut)
    transitionCounts.resize(options.maxCycleLength + 1);

  std::ostream &output = std::cout;
  if (options.outputType == options.probOut)
    output.precision(3);

  if (options.indexOfFirstNonOptionArgument == argc)
    processOneFile(std::cin, output);

  for (int i = options.indexOfFirstNonOptionArgument; i < argc; ++i) {
    izstream z;
    std::istream &input = openIn(argv[i], z);
    processOneFile(input, output);
  }

  if (options.outputType == options.countOut)
    writeCounts(output);

  return EXIT_SUCCESS;
}
catch( const std::bad_alloc& e ) {  // bad_alloc::what() may be unfriendly
  std::cerr << "tantan: out of memory\n";
  return EXIT_FAILURE;
}
catch (const std::exception &e) {
  std::cerr << "tantan: " << e.what() << '\n';
  return EXIT_FAILURE;
}
catch (int i) {
  return i;
}
