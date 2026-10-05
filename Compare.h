// =====================================================================
//  Compare.h
//  -------------------------------------------------------------------
//  LongPhase "compare" sub-command.
//
//  An independent C++ implementation; the metrics follow the definitions
//  used by whatshap compare so that results are directly comparable:
//      https://github.com/whatshap/whatshap
//  Phased VCFs are parsed with htslib and the per-block switch and
//  Hamming computations run in parallel.
// =====================================================================
#ifndef COMPARE_H
#define COMPARE_H

#include <string>

// Entry point dispatched from main.cpp
int CompareMain(int argc, char** argv, std::string version);

#endif
