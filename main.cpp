#include <fstream>
#include <iostream>
#include "Phasing.h"
#include "Haplotag.h"
#include "ModCall.h"
#include "Compare.h"

#ifdef HAS_ONNXRUNTIME
#include "GNN.h"
#endif

#define PROGRAM_BIN "longphase"
#define VERSION "2.0.2"

static std::string version = VERSION;

static const char *STRIDE_USAGE_MESSAGE =
"LongPhase Ver " VERSION "\n"
"Usage: " PROGRAM_BIN " <command> [options]\n"
"               phase      run phasing algorithm.\n"
"               haplotag   tag reads by haplotype.\n"
"               modcall    convert bam file to modification vcf file.\n"
"               compare    compare two phased VCF files.\n"
#ifdef HAS_ONNXRUNTIME
"               gnn        GNN-based post-hoc phasing correction.\n"
#endif
"\n";

int main(int argc, char** argv)
{
    if(argc <= 1)
    {
        std::cout << STRIDE_USAGE_MESSAGE;
        return 0;
    }

    std::string command(argv[1]);

    if(command=="phase")
    {
        PhasingMain(argc - 1, argv + 1, version);
    }
    else if(command=="haplotag")
    {
        HaplotagMain(argc - 1, argv + 1, version);
    }
    else if(command=="modcall")
    {
         ModCallMain(argc - 1, argv + 1, version);
    }
    else if(command=="compare")
    {
         CompareMain(argc - 1, argv + 1, version);
    }
#ifdef HAS_ONNXRUNTIME
    else if(command=="gnn")
    {
         GNNMain(argc - 1, argv + 1, version);
    }
#endif
    else{
        std::cout << STRIDE_USAGE_MESSAGE;
        return 0;
    }

    return 0;
}
