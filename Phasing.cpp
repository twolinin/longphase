#include "Phasing.h"
#include "PhasingProcess.h"
#include "Util.h"
#include <getopt.h>


#define SUBPROGRAM "phase"

static const char *CORRECT_USAGE_MESSAGE =
"Usage: "  " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"      --help                     display this help and exit.\n"
"      --dot                      each contig/chromosome will generate dot file. \n"
"      --ont, --pb                ont: Oxford Nanopore genomic reads.\n"
"                                 pb: PacBio HiFi/CCS genomic reads.\n"
"      --sv-file=NAME             input SV vcf file.\n"
"      -s, --snp-file=NAME        input SNP vcf file.\n"
"      -b, --bam-file=NAME        input bam file.\n"
"      -o, --out-prefix=NAME      prefix of phasing result.\n"
"      -r, --reference=NAME       reference fasta.\n"
"      -t, --threads=Num          number of thread. \n"
"      -d, --distance=Num         phasing two variant if distance less than threshold. default:300000\n"
"      -c, --crossBlock=Num       each block tries to connect with next N blocks. default:1\n"
"\n";


static const char* shortopts = "s:b:o:t:r:d:c:";

enum { OPT_HELP = 1 , DOT_FILE, SV_FILE, IS_ONT, IS_PB};

static const struct option longopts[] = { 
    { "help",              no_argument,        NULL, OPT_HELP },
    { "dot",               no_argument,        NULL, DOT_FILE },  
    { "ont",               no_argument,        NULL, IS_ONT }, 
    { "pb",                no_argument,        NULL, IS_PB }, 
    { "sv-file",           required_argument,  NULL, SV_FILE },     
    { "reference",         required_argument,  NULL, 'r' },
    { "snp-file",          required_argument,  NULL, 's' },
    { "bam-file",          required_argument,  NULL, 'b' },
    { "out-prefix",        required_argument,  NULL, 'o' },
    { "threads",           required_argument,  NULL, 't' },
    { "distance",          required_argument,  NULL, 'd' },
    { "crossBlock",        required_argument,  NULL, 'c' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static int numThreads = 1;
    static int distance = 300000;
    static int crossBlock = 1;
    static std::string snpFile="";
    static std::string svFile="";
    static std::string bamFile="";
    static std::string fastaFile="";
    static std::string resultPrefix="result";
    static bool generateDot=false;
    static bool isONT=false;
    static bool isPB=false;
}

void PhasingOptions(int argc, char** argv)
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
        case 'b': arg >> opt::bamFile; break;
        case 'o': arg >> opt::resultPrefix; break;
        case 'r': arg >> opt::fastaFile; break;  
        case 'd': arg >> opt::distance; break;  
        case 'c': arg >> opt::crossBlock; break;  
        case DOT_FILE: opt::generateDot=true; break;
        case IS_ONT: opt::isONT=true; break;
        case IS_PB: opt::isPB=true; break;
        case SV_FILE: arg >> opt::svFile; break; 
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
    
    if( opt::isONT == false && opt::isPB == false ){
        std::cerr << SUBPROGRAM ": missing arguments. --ont or --pb\n";
        die = true;
    }
    
    if( opt::isONT == true && opt::isPB == true ){
        std::cerr << SUBPROGRAM ": conflict arguments. --ont or --pb\n";
        die = true;
    }

    if( opt::snpFile != "")
    {
        std::ifstream openFile( opt::snpFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cout<< "File " << opt::snpFile << " not exist.\n\n";
            die = true;
        }
    }
    else{
        std::cerr << SUBPROGRAM ": missing SNP file.\n";
        die = true;
    }
    
    if( opt::bamFile != "")
    {
        std::ifstream openFile( opt::bamFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cout<< "File " << opt::bamFile << " not exist.\n\n";
            die = true;
        }
    }
    else{
        std::cerr << SUBPROGRAM ": missing bam file.\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << CORRECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }


    
}

int PhasingMain(int argc, char** argv)
{
    PhasingParameters ecParams;
    // set parameters
    PhasingOptions(argc, argv);
    // no file in command line

    ecParams.numThreads=opt::numThreads;
    ecParams.distance=opt::distance;
    ecParams.crossBlock=opt::crossBlock;
    ecParams.snpFile=opt::snpFile;
    ecParams.svFile=opt::svFile;
    ecParams.bamFile=opt::bamFile;
    ecParams.fastaFile=opt::fastaFile;
    ecParams.resultPrefix=opt::resultPrefix;
    ecParams.generateDot=opt::generateDot;
    ecParams.isONT=opt::isONT;
    ecParams.isPB=opt::isPB;
    
    PhasingProcess processor(ecParams);

    return 0;
}