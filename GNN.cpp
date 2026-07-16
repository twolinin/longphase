#include "GNN.h"
#include "GNNProcess.h"
#include <getopt.h>
#include <iostream>
#include <sstream>

#define SUBPROGRAM "gnn"

static const char *CORRECT_USAGE_MESSAGE =
"Usage: longphase " SUBPROGRAM " [OPTION]\n"
"      --help                          display this help and exit.\n\n"
"require arguments:\n"
"      -m, --model=NAME                input ONNX model file.\n"
"      -s, --vcf=NAME                  input phased VCF (SNV+Indel).\n"
"      --dot-prefix=NAME               DOT file prefix (from longphase phase --dot).\n"
"      -o, --output=NAME               output corrected VCF.\n"
"optional arguments:\n"
"      -r, --reference=NAME            reference FASTA (improves accuracy).\n"
"      --sv-vcf=NAME                   input phased SV VCF.\n"
"      --mod-vcf=NAME                  input phased methylation VCF.\n"
"      --output-sv-vcf=NAME            output corrected SV VCF.\n"
"      --output-mod-vcf=NAME           output corrected methylation VCF.\n"
"      -B, --break-threshold=FLOAT     unphase threshold. default:0.30\n"
"      --pe-threshold=FLOAT            PE trigger threshold. default:0.80\n"
"      --window=NUM                    window size in variants. default:20\n"
"      --respect-bridge                skip unphasing bridge vertices. default:false\n"
"      -t, --threads=NUM               number of threads. default:4\n";

static const char* shortopts = "m:s:o:r:B:t:";

enum {
    OPT_HELP = 1,
    OPT_DOT_PREFIX,
    OPT_SV_VCF,
    OPT_MOD_VCF,
    OPT_OUTPUT_SV_VCF,
    OPT_OUTPUT_MOD_VCF,
    OPT_PE_THRESHOLD,
    OPT_WINDOW,
    OPT_RESPECT_BRIDGE,
};

static const struct option longopts[] = {
    { "help",             no_argument,       NULL, OPT_HELP },
    { "model",            required_argument, NULL, 'm' },
    { "vcf",              required_argument, NULL, 's' },
    { "dot-prefix",       required_argument, NULL, OPT_DOT_PREFIX },
    { "output",           required_argument, NULL, 'o' },
    { "reference",        required_argument, NULL, 'r' },
    { "sv-vcf",           required_argument, NULL, OPT_SV_VCF },
    { "mod-vcf",          required_argument, NULL, OPT_MOD_VCF },
    { "output-sv-vcf",    required_argument, NULL, OPT_OUTPUT_SV_VCF },
    { "output-mod-vcf",   required_argument, NULL, OPT_OUTPUT_MOD_VCF },
    { "break-threshold",  required_argument, NULL, 'B' },
    { "pe-threshold",     required_argument, NULL, OPT_PE_THRESHOLD },
    { "window",           required_argument, NULL, OPT_WINDOW },
    { "respect-bridge",   no_argument,       NULL, OPT_RESPECT_BRIDGE },
    { "threads",          required_argument, NULL, 't' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static std::string modelFile = "";
    static std::string vcfFile = "";
    static std::string dotPrefix = "";
    static std::string outputVcf = "";
    static std::string referenceFile = "";
    static std::string svVcf = "";
    static std::string modVcf = "";
    static std::string outputSvVcf = "";
    static std::string outputModVcf = "";
    static float breakThreshold = 0.30f;
    static float peThreshold = 0.80f;
    static int window = 20;
    static int numThreads = 4;
    static bool respectBridge = false;
    static std::string command = "longphase ";
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
            case 'm':                arg >> opt::modelFile; break;
            case 's':                arg >> opt::vcfFile; break;
            case OPT_DOT_PREFIX:     arg >> opt::dotPrefix; break;
            case 'o':                arg >> opt::outputVcf; break;
            case 'r':                arg >> opt::referenceFile; break;
            case OPT_SV_VCF:        arg >> opt::svVcf; break;
            case OPT_MOD_VCF:       arg >> opt::modVcf; break;
            case OPT_OUTPUT_SV_VCF: arg >> opt::outputSvVcf; break;
            case OPT_OUTPUT_MOD_VCF:arg >> opt::outputModVcf; break;
            case 'B':                arg >> opt::breakThreshold; break;
            case OPT_PE_THRESHOLD:   arg >> opt::peThreshold; break;
            case OPT_WINDOW:         arg >> opt::window; break;
            case OPT_RESPECT_BRIDGE: opt::respectBridge = true; break;
            case 't':                arg >> opt::numThreads; break;
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

    if (opt::modelFile.empty()) {
        std::cerr << SUBPROGRAM ": missing model file (-m).\n";
        die = true;
    }
    if (opt::vcfFile.empty()) {
        std::cerr << SUBPROGRAM ": missing VCF file (-s).\n";
        die = true;
    }
    if (opt::outputVcf.empty()) {
        std::cerr << SUBPROGRAM ": missing output file (-o).\n";
        die = true;
    }
    if (opt::numThreads < 1) {
        std::cerr << SUBPROGRAM ": invalid threads. value: "
                  << opt::numThreads << "\n";
        die = true;
    }
    if (!opt::svVcf.empty() && opt::outputSvVcf.empty()) {
        std::cerr << SUBPROGRAM ": --output-sv-vcf required when --sv-vcf is given.\n";
        die = true;
    }
    if (!opt::modVcf.empty() && opt::outputModVcf.empty()) {
        std::cerr << SUBPROGRAM ": --output-mod-vcf required when --mod-vcf is given.\n";
        die = true;
    }

    // Auto-detect dot-prefix from VCF path
    if (opt::dotPrefix.empty()) {
        std::string vcf = opt::vcfFile;
        auto pos = vcf.rfind(".vcf");
        opt::dotPrefix = (pos != std::string::npos) ? vcf.substr(0, pos) : vcf;
        std::cerr << "[GNN] --dot-prefix not specified, using: "
                  << opt::dotPrefix << "\n";
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

    params.model_path      = opt::modelFile;
    params.vcf_path        = opt::vcfFile;
    params.dot_prefix      = opt::dotPrefix;
    params.reference_path  = opt::referenceFile;
    params.output_vcf      = opt::outputVcf;
    params.sv_vcf          = opt::svVcf;
    params.mod_vcf         = opt::modVcf;
    params.output_sv_vcf   = opt::outputSvVcf;
    params.output_mod_vcf  = opt::outputModVcf;
    params.break_threshold = opt::breakThreshold;
    params.pe_threshold    = opt::peThreshold;
    params.window          = opt::window;
    params.threads         = opt::numThreads;
    params.respect_bridge  = opt::respectBridge;

    std::cerr << "LongPhase Ver " << in_version << "\n\n"
              << "--- GNN Correction Parameters ---\n"
              << "Model           : " << params.model_path << "\n"
              << "Input VCF       : " << params.vcf_path << "\n"
              << "DOT prefix      : " << params.dot_prefix << "\n"
              << "Reference       : " << (params.reference_path.empty() ?
                                          "<not specified>" : params.reference_path) << "\n"
              << "Output VCF      : " << params.output_vcf << "\n"
              << "Break threshold : " << params.break_threshold << "\n"
              << "PE threshold    : " << params.pe_threshold << "\n"
              << "Window          : " << params.window << "\n"
              << "Respect bridge  : " << (params.respect_bridge ? "Yes" : "No") << "\n"
              << "Threads         : " << params.threads << "\n";
    if (!params.sv_vcf.empty())
        std::cerr << "SV VCF          : " << params.sv_vcf << "\n"
                  << "Output SV VCF   : " << params.output_sv_vcf << "\n";
    if (!params.mod_vcf.empty())
        std::cerr << "Methyl VCF      : " << params.mod_vcf << "\n"
                  << "Output Mod VCF  : " << params.output_mod_vcf << "\n";
    std::cerr << "\n";

    GNNModule module(params);
    return module.run();
}
