// Copyright 2010 Martin C. Frith

#include "mcf_tantan_options.hh"

#include "mcf_util.hh"

#include <unistd.h>

#include <cstdlib>  // EXIT_SUCCESS
#include <iostream>
#include <stdexcept>

namespace mcf {

typedef std::runtime_error Error;

static void badopt(char opt, const char *arg) {
  throw Error(std::string("bad option value: -") + opt + " " + arg);
}

static void writeAndQuit(const std::string &text) {
  std::cout << text;
  throw EXIT_SUCCESS;
}

static int myGetopt(int argc, char **argv, const char *optstring,
                    const std::string &help, const std::string &version) {
  if (optind < argc) {
    std::string nextarg = argv[optind];
    if (nextarg == "--help")    writeAndQuit(help);
    if (nextarg == "--version") writeAndQuit(version);
  }
  return getopt(argc, argv, optstring);
}

std::istream &operator>>(std::istream &s, TantanOptions::OutputType &x) {
  int i = 0;
  s >> i;
  if (i < 0 || i > 2)
    s.setstate(std::ios::failbit);
  if (s)
    x = static_cast<TantanOptions::OutputType>(i);
  return s;
}

TantanOptions::TantanOptions() :
    isProtein(false),
    maskSymbol(0),
    isPreserveLowercase(false),
    scoreMatrixFileName(0),
    repeatProb(0.005),
    repeatEndProb(0.05),
    maxCycleLength(-1),  // depends on isProtein
    repeatOffsetProbDecay(0.9),
    gapExistenceCost(0),
    gapExtensionCost(-1),  // means: no gaps
    minMaskProb(0.5),
    outputType(maskOut),
    indexOfFirstNonOptionArgument(-1) {}

void TantanOptions::fromArgs(int argc, char **argv) {
  std::string help = "\
Usage: tantan [options] fasta-sequence-file(s)\n\
Find simple repeats in sequences\n\
\n\
Options (default settings):\n\
 -p  interpret the sequences as proteins\n\
 -x  letter to use for masking, instead of lowercase\n\
 -c  preserve uppercase/lowercase in non-masked regions\n\
 -m  file for letter pair scores (+1/-1, but -p selects BLOSUM62)\n\
 -r  probability of a repeat starting per position ("
      + stringify(repeatProb) + ")\n\
 -e  probability of a repeat ending per position ("
      + stringify(repeatEndProb) + ")\n\
 -w  maximum tandem repeat period to consider (100, but -p selects 50)\n\
 -d  probability decay per period ("
      + stringify(repeatOffsetProbDecay) + ")\n\
 -a  gap existence cost ("
      + stringify(gapExistenceCost) + ")\n\
 -b  gap extension cost (infinite: no gaps)\n\
 -s  minimum repeat probability for masking ("
      + stringify(minMaskProb) + ")\n\
 -f  output type: 0=masked sequence, 1=repeat probabilities, 2=repeat counts ("
      + stringify(outputType) + ")\n\
 -h, --help  show help message, then exit\n\
 --version   show version information, then exit\n\
\n\
Report bugs to: tantan@cbrc.jp\n\
Home page: http://www.cbrc.jp/tantan/\n\
";

  std::string version = "tantan "
#include "version.hh"
      "\n";

  const char *optstring = "px:cm:r:e:w:d:a:b:s:f:h";

  int i;
  while ((i = myGetopt(argc, argv, optstring, help, version)) != -1) {
    char c = static_cast<char>(i);
    switch (c) {
      case 'p':
        isProtein = true;
        break;
      case 'x':
        unstringify(maskSymbol, optarg);
        break;
      case 'c':
        isPreserveLowercase = true;
        break;
      case 'm':
        scoreMatrixFileName = optarg;
        break;
      case 'r':
        unstringify(repeatProb, optarg);
        if (repeatProb < 0 || repeatProb >= 1)
          badopt(c, optarg);
        break;
      case 'e':
        unstringify(repeatEndProb, optarg);
        if (repeatEndProb < 0 || repeatEndProb > 1)
          badopt(c, optarg);
        break;
      case 'w':
        unstringify(maxCycleLength, optarg);
        if (maxCycleLength <= 0 )
          badopt(c, optarg);
        break;
      case 'd':
        unstringify(repeatOffsetProbDecay, optarg);
        if (repeatOffsetProbDecay <= 0 || repeatOffsetProbDecay > 1)
          badopt(c, optarg);
        break;
      case 'a':
        unstringify(gapExistenceCost, optarg);
        break;
      case 'b':
        unstringify(gapExtensionCost, optarg);
        if (gapExtensionCost <= 0)
          badopt(c, optarg);
        break;
      case 's':
        unstringify(minMaskProb, optarg);
        // don't bother checking for stupid values?
        break;
      case 'f':
        unstringify(outputType, optarg);
        break;
      case 'h':
        writeAndQuit(help);
      case '?':
        throw Error("bad option");
    }
  }

  if (gapExtensionCost > 0 && gapExistenceCost + gapExtensionCost <= 0)
    throw Error("gap existence + extension cost is too low");

  if (maxCycleLength < 0) maxCycleLength = (isProtein ? 50 : 100);

  indexOfFirstNonOptionArgument = optind;
}

}
