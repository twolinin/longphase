// =====================================================================
//  Compare.h
//  -------------------------------------------------------------------
//  LongPhase "compare" sub-command.
//
//  This module was adapted by Claude (Anthropic) from the WhatsHap
//  project, specifically:
//      https://github.com/whatshap/whatshap/blob/main/whatshap/cli/compare.py
//  It has been re-written in C++ with multi-thread support for the
//  switch / Hamming computations.
// =====================================================================
#ifndef COMPARE_H
#define COMPARE_H

#include <string>

// Entry point dispatched from main.cpp
int CompareMain(int argc, char** argv, std::string version);

#endif
