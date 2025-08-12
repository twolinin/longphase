#include "ModCall.h"
#include "ModCallProcess.h"
#include "Util.h"
#include <getopt.h>


#define SUBPROGRAM "modcall"

//
static const char *CORRECT_USAGE_MESSAGE =
"Usage: "  " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"      --help                        display this help and exit.\n"
"require arguments:\n"
"      -s, --snp-file=NAME           input SNP vcf file.\n"
"      -b, --bam-file=NAME           modified sorted bam file.\n"
"      -r, --reference=NAME          reference fasta.\n"

"optional arguments:\n"
"      -o, --out-prefix=NAME         prefix of phasing result. default: modcall_result\n"
"      -t, --threads=Num             number of thread. default:1\n"
"      --all                         output all coordinates where modifications have been detected.\n"

"phasing arguments:\n"
"      -m, --modThreshold=[0~1]      value extracted from MM tag and ML tag. \n"
"                                    above the threshold means modification occurred. default: 0.8\n"
"      -u, --unModThreshold=[0~1]    value extracted from MM tag and ML tag. \n"
"                                    above the threshold means no modification occurred. default: 0.2\n"
"      -e, --heterRatio=[0~1]        modified and unmodified scales. \n"
"                                    a higher threshold means that the two quantities need to be closer. default: 0.7\n"
"      -i, --noiseRatio=[0~1]        not being judged as modified and unmodified is noise.\n"
"                                    higher threshold means lower noise needs. default: 0.2\n"
"      -a, --connectAdjacent=Num     connect adjacent N METHs. default:6\n"
"      -c, --connectConfidence=[0~1] determine the confidence of phasing two ASMs.\n"
"                                    higher threshold requires greater consistency in the reads. default: 0.9\n"

"\n";

static const char* shortopts = "s:o:t:r:b:m:u:e:i:";

enum { OPT_HELP = 1, ALL_MOD };

static const struct option longopts[] = { 
    { "help",              no_argument,        NULL, OPT_HELP }, 
    { "all",               no_argument,        NULL, ALL_MOD },
    { "reference",         required_argument,  NULL, 'r' },
    { "snp-file",          required_argument,  NULL, 's' },
    { "out-prefix",        required_argument,  NULL, 'o' },
    { "threads",           required_argument,  NULL, 't' },
	{ "methylbamfile",     required_argument,  NULL, 'b' },
    { "modThreshold",      required_argument,  NULL, 'm' },
    { "unModThreshold",    required_argument,  NULL, 'u' },
    { "heterRatio",        required_argument,  NULL, 'e' },
    { "noiseRatio",        required_argument,  NULL, 'i' },
    { "connectAdjacent",   required_argument,  NULL, 'a' },
    { "connectConfidence", required_argument,  NULL, 'c' },
    { "iterCount",         required_argument,  NULL, 'k' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static int numThreads = 1;
    static std::string fastaFile = "";
    static std::string snpFile="";
    static std::string resultPrefix = "modcall_result";
    static std::vector<std::string> bamFileVec;
    static float modThreshold = 0.8;
    static float unModThreshold = 0.2;
    static float heterRatio = 0.6;
    static float noiseRatio = 0.2;
    static int connectAdjacent = 20;
    static float connectConfidence = 0.9;
    static int iterCount = 2;
    static std::string command;
    
    bool outputAllMod = false;
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
        case 's': arg >> opt::snpFile; break;
        case 't': arg >> opt::numThreads; break;
        case 'o': arg >> opt::resultPrefix; break;
        case 'r': arg >> opt::fastaFile; break;  
        case 'm': arg >> opt::modThreshold; break;
        case 'u': arg >> opt::unModThreshold; break;
        case 'e': arg >> opt::heterRatio; break;
        case 'i': arg >> opt::noiseRatio; break;
        case 'a': arg >> opt::connectAdjacent; break;
        case 'c': arg >> opt::connectConfidence; break;
        case 'k': arg >> opt::iterCount; break;
        case 'b': {
            std::string bamFile;
            arg >> bamFile;
            opt::bamFileVec.push_back(bamFile); break;
        }
        case ALL_MOD: opt::outputAllMod=true; break;
        case OPT_HELP:
            std::cout << CORRECT_USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        }
    }
    
    for(int i = 0; i < argc; ++i){
        opt::command.append(argv[i]);
        opt::command.append(" ");
    }

    if (argc - optind < 0 )
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }

    if( opt::snpFile != "")
    {
        std::ifstream openFile( opt::snpFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cerr<< "File " << opt::snpFile << " not exist.\n\n";
            die = true;
        }
    }
    
	if( opt::bamFileVec.size() != 0 )
    {
        
        for (auto bamFile: opt::bamFileVec ){
            std::ifstream openFile( bamFile.c_str() );
            if( !openFile.is_open() )
            {
                std::cout<< "File " << bamFile << " not exist.\n\n";
                die = true;
            }
        }
    }
    else{
        std::cerr << SUBPROGRAM ": missing methylBamFile file.\n";
        die = true;
    }
    
    if( opt::fastaFile != "")
    {
        std::ifstream openFile( opt::fastaFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cerr<< "File " << opt::fastaFile << " not exist.\n\n";
            die = true;
        }
    }
    else{
        std::cerr << SUBPROGRAM ": missing reference.\n";
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

int ModCallMain(int argc, char** argv, std::string in_version)
{
    ModCallParameters ecParams;
    // set parameters
    ModCallOptions(argc, argv);

    ecParams.numThreads=opt::numThreads;
    ecParams.fastaFile=opt::fastaFile;
    ecParams.snpFile=opt::snpFile;
    ecParams.resultPrefix=opt::resultPrefix;
    ecParams.bamFileVec=opt::bamFileVec;
    ecParams.modThreshold=opt::modThreshold;
    ecParams.unModThreshold=opt::unModThreshold;
    ecParams.heterRatio=opt::heterRatio;
    ecParams.noiseRatio=opt::noiseRatio;
    ecParams.connectAdjacent=opt::connectAdjacent;
    ecParams.connectConfidence=opt::connectConfidence;
    ecParams.iterCount=opt::iterCount;
    ecParams.version=in_version;
    ecParams.command=opt::command;
    
    ecParams.outputAllMod=opt::outputAllMod;
    
    ModCallProcess processor(ecParams);

    return 0;
}
