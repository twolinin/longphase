// =====================================================================
//  Compare.cpp
//  -------------------------------------------------------------------
//  LongPhase "compare" sub-command – command-line front-end.
//
//  Adapted by Claude (Anthropic) from WhatsHap's compare module:
//      https://github.com/whatshap/whatshap/blob/main/whatshap/cli/compare.py
//
//  This file is only responsible for parsing the CLI arguments and
//  invoking CompareProcess.  All of the real work (VCF parsing, block
//  intersection, switch / Hamming computation, threading, output) is
//  in CompareProcess.cpp.
// =====================================================================

#include "Compare.h"
#include "CompareProcess.h"

#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// ------------------------------------------------------------------ //
// Usage text
// ------------------------------------------------------------------ //
static const char *COMPARE_USAGE =
"Usage:  compare [OPTION] ... TRUTH.vcf  QUERY1.vcf [QUERY2.vcf ...]\n"
"        Compare two or more phased VCF files (first VCF is treated as truth).\n"
"\n"
"  -h, --help                     display this help and exit.\n"
"\n"
"require arguments:\n"
"      At least two phased VCF/BCF files.  The first one is treated\n"
"      as the ground-truth phasing.\n"
"\n"
"optional arguments:\n"
"  -o, --out-prefix=NAME          prefix of output TSV file.\n"
"                                 default: compare_result\n"
"  -t, --threads=NUM              number of threads used for the parallel\n"
"                                 switch / Hamming computation. default: 1\n"
"  -n, --names=N1,N2,...          comma-separated list of data-set names\n"
"                                 used in the report (same order as VCFs).\n"
"  -s, --sample=SAMPLE            name of the sample to process.\n"
"                                 default: first sample found in VCF\n"
"      --ignore-sample-name       for single-sample VCFs, ignore sample\n"
"                                 name and assume all samples are the same.\n"
"      --only-snvs                only process SNVs, ignore all other\n"
"                                 variants.\n"
"      --sw-bed=FILE              write per-switch-error BED records to\n"
"                                 FILE (useful for debugging).  Columns:\n"
"                                 chrom, start(0-based), end(exclusive),\n"
"                                 type, ps_truth, ps_query, dataset_pair.\n"
"\n";

// ------------------------------------------------------------------ //
// Helper: split a comma-separated string
// ------------------------------------------------------------------ //
static std::vector<std::string> splitComma(const std::string& s)
{
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) out.push_back(item);
    return out;
}

// ------------------------------------------------------------------ //
// Entry point – called from main.cpp
// ------------------------------------------------------------------ //
int CompareMain(int argc, char** argv, std::string version)
{
    CompareParameters params;
    params.version      = version;
    params.outPrefix    = "compare_result";
    params.numThreads   = 1;
    params.sample       = "";
    params.ignoreSampleName = false;
    params.onlySnvs     = false;

    // Remember the full command-line for the ##commandline= header
    {
        std::ostringstream cmd;
        for (int i = 0; i < argc; ++i) {
            if (i) cmd << ' ';
            cmd << argv[i];
        }
        params.command = cmd.str();
    }

    static struct option long_opts[] = {
        {"help",               no_argument,       0, 'h'},
        {"out-prefix",         required_argument, 0, 'o'},
        {"threads",            required_argument, 0, 't'},
        {"names",              required_argument, 0, 'n'},
        {"sample",             required_argument, 0, 's'},
        {"ignore-sample-name", no_argument,       0,  1 },
        {"only-snvs",          no_argument,       0,  2 },
        {"sw-bed",             required_argument, 0,  3 },
        {0, 0, 0, 0}
    };

    int opt, idx;
    std::string namesArg;
    while ((opt = getopt_long(argc, argv, "ho:t:n:s:", long_opts, &idx)) != -1) {
        switch (opt) {
            case 'h': std::cout << COMPARE_USAGE; return 0;
            case 'o': params.outPrefix = optarg;                 break;
            case 't': params.numThreads = std::atoi(optarg);     break;
            case 'n': namesArg = optarg;                         break;
            case 's': params.sample = optarg;                    break;
            case  1 : params.ignoreSampleName = true;            break;
            case  2 : params.onlySnvs = true;                    break;
            case  3 : params.swBedFile = optarg;                 break;
            default : std::cout << COMPARE_USAGE; return 1;
        }
    }

    // Positional = VCF files
    while (optind < argc) params.vcfFiles.emplace_back(argv[optind++]);

    if (params.vcfFiles.size() < 2) {
        std::cerr << "[compare] ERROR: at least two VCF files are required.\n\n";
        std::cout << COMPARE_USAGE;
        return 1;
    }

    if (params.numThreads < 1) params.numThreads = 1;
    unsigned hw = std::thread::hardware_concurrency();
    if (hw && (unsigned)params.numThreads > hw * 4) {
        std::cerr << "[compare] Warning: requested threads (" << params.numThreads
                  << ") much larger than hardware_concurrency (" << hw << ").\n";
    }

    // Dataset names
    if (!namesArg.empty()) {
        params.datasetNames = splitComma(namesArg);
        if (params.datasetNames.size() != params.vcfFiles.size()) {
            std::cerr << "[compare] ERROR: number of --names entries ("
                      << params.datasetNames.size()
                      << ") does not match number of VCFs ("
                      << params.vcfFiles.size() << ")\n";
            return 1;
        }
    } else {
        for (size_t i = 0; i < params.vcfFiles.size(); ++i)
            params.datasetNames.emplace_back("file" + std::to_string(i));
    }

    // Print parameter summary (same style as PhasingProcess.cpp)
    std::cerr << "LongPhase Ver " << params.version << "\n\n";
    std::cerr << "--- compare Parameter ---\n";
    std::cerr << "Output Prefix      : " << params.outPrefix  << "\n";
    std::cerr << "Number of Threads  : " << params.numThreads << "\n";
    std::cerr << "Only SNVs          : " << (params.onlySnvs ? "True" : "False") << "\n";
    std::cerr << "Sample             : " << (params.sample.empty() ? "<first>" : params.sample) << "\n";
    std::cerr << "SW BED             : " << (params.swBedFile.empty() ? "<disabled>" : params.swBedFile) << "\n";
    std::cerr << "VCF files          :\n";
    for (size_t i = 0; i < params.vcfFiles.size(); ++i)
        std::cerr << "   [" << params.datasetNames[i] << "] " << params.vcfFiles[i] << "\n";
    std::cerr << "\n";

    // Dispatch to the real work
    CompareProcess cp(params);
    return cp.run();
}
