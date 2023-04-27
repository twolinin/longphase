#include "ModCall.h"
#include "ModCallProcess.h"
#include "Util.h"
#include <getopt.h>


#define SUBPROGRAM "modcall"

//
static const char *CORRECT_USAGE_MESSAGE =
"Usage: "  " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"      --help                        display this help and exit.\n"
"      -o, --out-prefix=NAME         prefix of phasing result.\n"
"      -r, --reference=NAME          reference fasta. default: 1\n"
"      -t, --threads=Num             number of thread. \n"
"      -b, --bam-file=NAME           modified sorted bam file.\n"
"      -m, --modThreshold=[0~1]      value extracted from MM tag and ML tag. \n"
"                                    above the threshold means modification occurred. default: 0.8\n"
"      -u, --unModThreshold=[0~1]    value extracted from MM tag and ML tag. \n"
"                                    above the threshold means no modification occurred. default: 0.2\n"
"      -e, --heterRatio=[0~1]        modified and unmodified scales. \n"
"                                    a higher threshold means that the two quantities need to be closer. default: 0.5\n"
"      -i, --noiseRatio=[0~1]        not being judged as modified and unmodified is noise.\n"
"                                    higher threshold means lower noise needs. default: 0.3\n"
"\n";

static const char* shortopts = "o:t:r:b:m:u:e:i:";

enum { OPT_HELP = 1 };

static const struct option longopts[] = { 
    { "help",              no_argument,        NULL, OPT_HELP }, 
    { "reference",         required_argument,  NULL, 'r' },
    { "out-prefix",        required_argument,  NULL, 'o' },
    { "threads",           required_argument,  NULL, 't' },
	{ "methylbamfile",     required_argument,  NULL, 'b' },
    { "modThreshold",      required_argument,  NULL, 'm' },
    { "unModThreshold",    required_argument,  NULL, 'u' },
    { "heterRatio",        required_argument,  NULL, 'e' },
    { "noiseRatio",        required_argument,  NULL, 'i' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static int numThreads = 1;
    static std::string fastaFile = "";
    static std::string resultPrefix = "result";
	static std::string methylBamFile = "";
    static float modThreshold = 0.8;
    static float unModThreshold = 0.2;
    static float heterRatio = 0.5;
    static float noiseRatio = 0.3;
}

void ModCallOptions(int argc, char** argv)
{
    optind=1;    //reset getopt
    
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
        case 't': arg >> opt::numThreads; break;
        case 'o': arg >> opt::resultPrefix; break;
        case 'r': arg >> opt::fastaFile; break;  
		case 'b': arg >> opt::methylBamFile; break;		
        case 'm': arg >> opt::modThreshold; break;
        case 'u': arg >> opt::unModThreshold; break;
        case 'e': arg >> opt::heterRatio; break;
        case 'i': arg >> opt::noiseRatio; break;
        case OPT_HELP:
            std::cout << CORRECT_USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 0 )
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }

	if( opt::methylBamFile != "")
    {
        std::ifstream openFile( opt::methylBamFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cout<< "File " << opt::methylBamFile << " not exist.\n\n";
            die = true;
        }
    }
    else{
        std::cerr << SUBPROGRAM ": missing methylBamFile file.\n";
        die = true;
    }
    
    if( opt::modThreshold < opt::unModThreshold ){
        std::cerr << "error: modThreshold is lower than unModThreshold. Please Check -m and -u\n";
        std::cerr << "modThreshold: " << opt::modThreshold << "\n";
        std::cerr << "unModThreshold: " << opt::unModThreshold << "\n";
        die = true;
    }
    
    if (die)
    {
        std::cout << "\n" << CORRECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

}

int ModCallMain(int argc, char** argv)
{
    ModCallParameters ecParams;
    // set parameters
    ModCallOptions(argc, argv);

    ecParams.numThreads=opt::numThreads;
    ecParams.fastaFile=opt::fastaFile;
    ecParams.resultPrefix=opt::resultPrefix;
	ecParams.methylBamFile=opt::methylBamFile;
    ecParams.modThreshold=opt::modThreshold;
    ecParams.unModThreshold=opt::unModThreshold;
    ecParams.heterRatio=opt::heterRatio;
    ecParams.noiseRatio=opt::noiseRatio;
    
    ModCallProcess processor(ecParams);

    return 0;
}