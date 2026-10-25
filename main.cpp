#include <fstream>
#include <iostream>
#include "Phasing.h"
#include "Haplotag.h"
#include "ModCall.h"
#include "Compare.h"


#define PROGRAM_BIN "longphase"
#define VERSION "2.1.2"



static std::string version = VERSION;

// jemalloc: purge freed pages from a background thread. Without it, memory freed by a finished
// chromosome stays resident in that thread's arena, so the peak is the sum of all chromosomes
// rather than of those in progress. MALLOC_CONF in the environment still overrides this.
extern "C" {
const char *malloc_conf = "background_thread:true";
}

static const char *STRIDE_USAGE_MESSAGE =
"LongPhase Ver " VERSION "\n"
"Usage: " PROGRAM_BIN " <command> [options]\n"  
"               phase      run phasing algorithm.\n"
"               haplotag   tag reads by haplotype.\n"
"               modcall    convert bam file to modification vcf file.\n"
"               compare    compare two phased VCF files.\n"

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
    else{
        std::cout << STRIDE_USAGE_MESSAGE;
        return 0;
    }

    return 0;
}
