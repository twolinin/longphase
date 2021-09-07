#include <fstream>
#include <iostream>
#include "Phasing.h"

#define PROGRAM_BIN "main"

static const char *STRIDE_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " <command> [options]\n"  
"               phase      run phasing algorithm.\n"
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
        PhasingMain(argc - 1, argv + 1);
    }
    else{
        std::cout << STRIDE_USAGE_MESSAGE;
        return 0;
    }

    return 0;
}
