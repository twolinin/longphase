#include <fstream>
#include <iostream>
#include "Phasing.h"
#include "Haplotag.h"
#include "ModCall.h"


#define PROGRAM_BIN "main"
#define VERSION "1.7.2dev"

static std::string version = VERSION;

static const char *STRIDE_USAGE_MESSAGE =
"Version: " VERSION " \n"
"Usage: " PROGRAM_BIN " <command> [options]\n"  
"               phase      run phasing algorithm.\n"
"               haplotag   tag reads by haplotype.\n"
"               modcall    convert bam file to modification vcf file.\n"

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
        HaplotagMain(argc - 1, argv + 1);
    }
    else if(command=="modcall")
    {
         ModCallMain(argc - 1, argv + 1, version);
    }
    else{
        std::cout << STRIDE_USAGE_MESSAGE;
        return 0;
    }

    return 0;
}
