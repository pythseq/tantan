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
  int gapExistenceCost;
  int gapExtensionCost;
  double minMaskProb;
  enum OutputType { maskOut, probOut, countOut, bedOut } outputType;

  int indexOfFirstNonOptionArgument;
};

}

#endif
