// Copyright 2010 Martin C. Frith

// Mask simple regions (low complexity & short-period tandem repeats)
// in biological sequences.

#include "mcf_alphabet.hh"
#include "mcf_fasta_sequence.hh"
#include "mcf_score_matrix.hh"
#include "mcf_score_matrix_probs.hh"
#include "mcf_tantan_options.hh"
#include "mcf_util.hh"
#include "tantan.hh"

#include <algorithm>  // fill_n
#include <cassert>
#include <cmath>
#include <cstring>  // strchr
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream>
#include <iostream>

#define BEG(v) ((v).empty() ? 0 : &(v).front())
#define END(v) ((v).empty() ? 0 : &(v).back() + 1)

typedef std::runtime_error Error;

using namespace mcf;

// move this function to a reusable file?
bool isDubiousDna(const uchar *beg, const uchar *end) {
  int badLetters = 0;
  if (end - beg > 100) end = beg + 100;  // just check the first 100 letters
  while (beg < end) {
    if (!std::strchr("AaCcGgTtNnUu", *beg)) {
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

  if (options.scoreMatrixFileName)
    unfilify(scoreMatrix, options.scoreMatrixFileName);
  else if (options.isProtein)
    unstringify(scoreMatrix, ScoreMatrix::blosum62);
  else
    scoreMatrix.initMatchMismatch(1, 1, "ACGTU");  // allow for RNA

  scoreMatrix.makeFastMatrix(fastMatrixPointers, scoreMatrixSize,
                             alphabet.lettersToNumbers,
                             scoreMatrix.minScore(), false);

  ScoreMatrixProbs smp(fastMatrixPointers, alphabet.size);
  if (smp.isBad())
    throw Error("can't calculate probabilities for this score matrix");

  for (int i = 0; i < scoreMatrixSize; ++i)
    for (int j = 0; j < scoreMatrixSize; ++j)
      probMatrix[i][j] = std::exp(smp.lambda() * fastMatrix[i][j]);

  if (options.gapExtensionCost > 0) {
    int firstGapCost = options.gapExistenceCost + options.gapExtensionCost;
    firstGapProb = std::exp(-smp.lambda() * firstGapCost);
    otherGapProb = std::exp(-smp.lambda() * options.gapExtensionCost);

    // the gap existence cost includes the gap ending cost:
    firstGapProb /= (1 - otherGapProb);

    // XXX check if firstGapProb is too high
  }

  //std::cerr << "lambda: " << smp.lambda() << "\n";
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
  } else if(options.outputType == options.probOut) {
    std::vector<float> probabilities(end - beg);
    float *probBeg = BEG(probabilities);
    float *probEnd = END(probabilities);
    tantan::getProbabilities(beg, end, options.maxCycleLength,
                             probMatrixPointers,
                             options.repeatProb, options.repeatEndProb,
                             options.repeatOffsetProbDecay,
                             firstGapProb, otherGapProb, probBeg);
    output << '>' << f.title << '\n';
    for (float *i = probBeg; i < probEnd; ++i)
      output << *i << '\n';
  } else {
    tantan::countTransitions(beg, end, options.maxCycleLength,
                             probMatrixPointers,
                             options.repeatProb, options.repeatEndProb,
                             options.repeatOffsetProbDecay,
                             firstGapProb, otherGapProb,
                             BEG(transitionCounts));
    double sequenceLength = static_cast<double>(end - beg);
    transitionTotal += sequenceLength + 1;
  }
}

void processOneFile(std::istream &input, std::ostream &output) {
  bool isFirstSequence = true;

  // This code strives to minimize memory usage.  The sequence-reading
  // operation does not overwrite the sequence until it finishes
  // reading successfully.  So, we don't want to overwrite a large,
  // old sequence.  Hence, we make a brand-new FastaSequence each time
  // through the loop.
  while (true) {
    FastaSequence f;
    if (!(input >> f)) break;
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
    std::ifstream ifs;
    std::istream &input = openIn(argv[i], ifs);
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
