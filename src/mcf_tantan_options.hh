// Copyright 2010 Martin C. Frith

#ifndef MCF_TANTAN_OPTIONS_HH
#define MCF_TANTAN_OPTIONS_HH

namespace mcf {

struct TantanOptions {
  TantanOptions();

  void fromArgs(int argc, char** argv);

  bool isProtein;
  char maskSymbol;
  bool isPreserveLowercase;
  const char *scoreMatrixFileName;
  double repeatProb;
  double repeatEndProb;
  int maxCycleLength;
  double repeatOffsetProbDecay;
  int matchScore;
  int mismatchCost;
  int gapExistenceCost;
  int gapExtensionCost;
  double minMaskProb;
  double minCopyNumber;
  enum OutputType { maskOut, probOut, countOut, bedOut, repOut } outputType;

  int indexOfFirstNonOptionArgument;
};

}

#endif
