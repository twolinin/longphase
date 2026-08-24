#include "GNN.h"
#include "GNNProcess.h"
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <glob.h>

#define SUBPROGRAM "gnn"

static const char *CORRECT_USAGE_MESSAGE =
"Usage: longphase " SUBPROGRAM " [OPTION]\n"
"      --help                          display this help and exit.\n\n"
"require arguments:\n"
"      -s, --snp-file=NAME             input phased SNP/SNV vcf file (from longphase phase).\n"
"      -o, --out-prefix=NAME           prefix of corrected result. default:result\n"
"optional arguments:\n"
"      -r, --reference=NAME            reference fasta. improves accuracy.\n"
"      --sv-file=NAME                  input phased SV vcf file.\n"
"      --mod-file=NAME                 input phased modified vcf file.\n"
"      --dot-prefix=NAME               DOT file prefix. default: same as --snp-file without .vcf\n"
"      -B, --break-threshold=[0~1]     unphase a variant when GNN error probability exceeds\n"
"                                      this value. default:0.30\n"
"      --pe-threshold=[0~1]            phasing entropy threshold to trigger GNN. default:0.80\n"
"      --window=NUM                    window size in variants. default:20\n"
"      --respect-bridge                do not unphase bridge vertices. default:false\n"
"      --no-split-blocks               do not split PS blocks that become disconnected\n"
"                                      after a bridge variant is unphased. default:false\n"
"      -t, --threads=NUM               number of thread. default:4\n"
"\n"
"Output files (based on -o prefix):\n"
"      <prefix>.vcf                    corrected SNP vcf\n"
"      <prefix>_SV.vcf                 corrected SV vcf   (only if --sv-file given)\n"
"      <prefix>_mod.vcf                corrected mod vcf  (only if --mod-file given)\n";

static const char* shortopts = "s:o:r:B:t:";

enum {
    OPT_HELP = 1,
    OPT_DOT_PREFIX,
    OPT_SV_FILE,
    OPT_MOD_FILE,
    OPT_PE_THRESHOLD,
    OPT_WINDOW,
    OPT_RESPECT_BRIDGE,
    OPT_NO_SPLIT_BLOCKS,
};

static const struct option longopts[] = {
    { "help",             no_argument,       NULL, OPT_HELP },
    { "snp-file",         required_argument, NULL, 's' },
    { "out-prefix",       required_argument, NULL, 'o' },
    { "reference",        required_argument, NULL, 'r' },
    { "sv-file",          required_argument, NULL, OPT_SV_FILE },
    { "mod-file",         required_argument, NULL, OPT_MOD_FILE },
    { "dot-prefix",       required_argument, NULL, OPT_DOT_PREFIX },
    { "break-threshold",  required_argument, NULL, 'B' },
    { "pe-threshold",     required_argument, NULL, OPT_PE_THRESHOLD },
    { "window",           required_argument, NULL, OPT_WINDOW },
    { "respect-bridge",   no_argument,       NULL, OPT_RESPECT_BRIDGE },
    { "no-split-blocks",  no_argument,       NULL, OPT_NO_SPLIT_BLOCKS },
    { "threads",          required_argument, NULL, 't' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static std::string snpFile = "";
    static std::string resultPrefix = "result";
    static std::string referenceFile = "";
    static std::string svFile = "";
    static std::string modFile = "";
    static std::string dotPrefix = "";
    static float breakThreshold = 0.30f;
    static float peThreshold = 0.80f;
    static int window = 20;
    static int numThreads = 4;
    static bool respectBridge = false;
    static bool splitBlocks = true;
    static std::string command = "longphase ";
}

// Count files matching "<prefix>.*.dot"
static int countDotFiles(const std::string& prefix)
{
    glob_t g;
    std::string pattern = prefix + ".*.dot";
    int n = 0;
    if (glob(pattern.c_str(), 0, NULL, &g) == 0) {
        n = (int)g.gl_pathc;
    }
    globfree(&g);
    return n;
}

// Strip a trailing .vcf / .vcf.gz / .VCF from a path
static std::string stripVcfSuffix(const std::string& path)
{
    static const char* sfx[] = { ".vcf.gz", ".vcf", ".VCF", ".bcf" };
    for (const char* s : sfx) {
        size_t sl = strlen(s);
        if (path.size() > sl && path.compare(path.size() - sl, sl, s) == 0)
            return path.substr(0, path.size() - sl);
    }
    return path;
}

void GNNOptions(int argc, char** argv)
{
    optind = 1;

    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 's': arg >> opt::snpFile; break;
            case 'o': arg >> opt::resultPrefix; break;
            case 'r': arg >> opt::referenceFile; break;
            case 'B': arg >> opt::breakThreshold; break;
            case 't': arg >> opt::numThreads; break;
            case OPT_SV_FILE:         arg >> opt::svFile; break;
            case OPT_MOD_FILE:        arg >> opt::modFile; break;
            case OPT_DOT_PREFIX:      arg >> opt::dotPrefix; break;
            case OPT_PE_THRESHOLD:    arg >> opt::peThreshold; break;
            case OPT_WINDOW:          arg >> opt::window; break;
            case OPT_RESPECT_BRIDGE:  opt::respectBridge = true; break;
            case OPT_NO_SPLIT_BLOCKS: opt::splitBlocks = false; break;
            case OPT_HELP:
                std::cout << CORRECT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            default: die = true; break;
        }
    }

    for (int i = 0; i < argc; ++i) {
        opt::command.append(argv[i]);
        opt::command.append(" ");
    }

    if (opt::snpFile != "") {
        std::ifstream openFile(opt::snpFile.c_str());
        if (!openFile.is_open()) {
            std::cerr << "File " << opt::snpFile << " not exist.\n\n";
            die = true;
        }
    }
    else {
        std::cerr << SUBPROGRAM ": missing SNP file.\n";
        die = true;
    }

    if (opt::svFile != "") {
        std::ifstream openFile(opt::svFile.c_str());
        if (!openFile.is_open()) {
            std::cerr << "File " << opt::svFile << " not exist.\n\n";
            die = true;
        }
    }

    if (opt::modFile != "") {
        std::ifstream openFile(opt::modFile.c_str());
        if (!openFile.is_open()) {
            std::cerr << "File " << opt::modFile << " not exist.\n\n";
            die = true;
        }
    }

    if (opt::numThreads < 1) {
        std::cerr << SUBPROGRAM " invalid threads. value: "
                  << opt::numThreads
                  << "\nplease check -t, --threads=Num\n";
        die = true;
    }

    if (opt::breakThreshold < 0 || opt::breakThreshold > 1) {
        std::cerr << SUBPROGRAM " invalid break threshold. value: "
                  << opt::breakThreshold
                  << "\nthis value need: 0~1, please check -B, --break-threshold=[0~1]\n";
        die = true;
    }

    if (opt::peThreshold < 0 || opt::peThreshold > 1) {
        std::cerr << SUBPROGRAM " invalid PE threshold. value: "
                  << opt::peThreshold
                  << "\nthis value need: 0~1, please check --pe-threshold=[0~1]\n";
        die = true;
    }

    if (opt::window < 1) {
        std::cerr << SUBPROGRAM " invalid window. value: "
                  << opt::window
                  << "\nplease check --window=NUM\n";
        die = true;
    }

    // ── Resolve DOT prefix ───────────────────────────────
    // Default: same path as --snp-file with the .vcf suffix removed.
    if (!die) {
        if (opt::dotPrefix == "") {
            std::string guess = stripVcfSuffix(opt::snpFile);
            if (countDotFiles(guess) > 0) {
                opt::dotPrefix = guess;
            }
            else {
                std::cerr << SUBPROGRAM ": no DOT file found at \""
                          << guess << ".*.dot\".\n"
                          << "DOT files are produced by 'longphase phase --dot'.\n"
                          << "If they are stored elsewhere, specify the prefix with"
                          << " --dot-prefix=NAME\n";
                die = true;
            }
        }
        else if (countDotFiles(opt::dotPrefix) == 0) {
            std::cerr << SUBPROGRAM ": no DOT file found at \""
                      << opt::dotPrefix << ".*.dot\".\n"
                      << "please check --dot-prefix=NAME\n";
            die = true;
        }
    }

    if (die) {
        std::cerr << "\n" << CORRECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int GNNMain(int argc, char** argv, std::string in_version)
{
    GNNModule::Params params;

    GNNOptions(argc, argv);

    params.vcf_path        = opt::snpFile;
    params.dot_prefix      = opt::dotPrefix;
    params.reference_path  = opt::referenceFile;
    params.sv_vcf          = opt::svFile;
    params.mod_vcf         = opt::modFile;
    params.break_threshold = opt::breakThreshold;
    params.pe_threshold    = opt::peThreshold;
    params.window          = opt::window;
    params.threads         = opt::numThreads;
    params.respect_bridge  = opt::respectBridge;
    params.split_blocks    = opt::splitBlocks;

    // Derive output paths from the prefix; secondary outputs are only
    // produced when the matching input was supplied.
    params.output_vcf = opt::resultPrefix + ".vcf";
    if (!params.sv_vcf.empty())
        params.output_sv_vcf = opt::resultPrefix + "_SV.vcf";
    if (!params.mod_vcf.empty())
        params.output_mod_vcf = opt::resultPrefix + "_mod.vcf";

    std::cerr << "LongPhase Ver " << in_version << "\n\n"
              << "--- GNN Correction Parameters ---\n"
              << "Model           : built-in\n"
              << "SNP file        : " << params.vcf_path << "\n";
    if (!params.sv_vcf.empty())
        std::cerr << "SV file         : " << params.sv_vcf << "\n";
    if (!params.mod_vcf.empty())
        std::cerr << "Mod file        : " << params.mod_vcf << "\n";
    std::cerr << "DOT prefix      : " << params.dot_prefix << "\n"
              << "Reference       : " << (params.reference_path.empty() ?
                                          "<not specified>" : params.reference_path) << "\n"
              << "Output prefix   : " << opt::resultPrefix << "\n"
              << "Break threshold : " << params.break_threshold << "\n"
              << "PE threshold    : " << params.pe_threshold << "\n"
              << "Window          : " << params.window << "\n"
              << "Respect bridge  : " << (params.respect_bridge ? "Yes" : "No") << "\n"
              << "Split blocks    : " << (params.split_blocks ? "Yes" : "No") << "\n"
              << "Threads         : " << params.threads << "\n\n";

    GNNModule module(params);
    return module.run();
}
