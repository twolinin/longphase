#include "Phasing.h"
#include "PhasingProcess.h"
#include "Util.h"
#include <getopt.h>


#define SUBPROGRAM "phase"

static const char *CORRECT_USAGE_MESSAGE =
"Usage: "  " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"      --help                          display this help and exit.\n\n"
"require arguments:\n"
"      -s, --snp-file=NAME             input SNP vcf file.\n"
"      -b, --bam-file=NAME             input bam file.\n"
"      --ont, --pb                     ont: Oxford Nanopore genomic reads.\n"
"                                      pb: PacBio HiFi/CCS genomic reads.\n"
"optional arguments:\n"
"      -r, --reference=NAME            reference fasta.\n"
"      --sv-file=NAME                  input SV vcf file.\n"
"      -t, --threads=Num               number of thread. default:1\n"
"      -o, --out-prefix=NAME           prefix of phasing result.\n"
"      --dot                           each contig/chromosome will generate dot file. \n"
"      -1, --readsThreshold=[0~1]      give up SNP-SNP phasing pair if the number of reads of the two combinations are similar. default:0.5\n"
"      -2, --qualityThreshold=[0~1]    give up phasing pair if the quality of the two combinations are similar. default:0.7\n"
"      -3, --blockReadThreshold=[0~1]  give up phasing pair if the number of block's reads of the two combinations are similar. default:1.0\n"
"      -4, --svReadsThreshold=[0~1]    give up SNP-SV phasing pair if the number of reads of the two combinations are similar. default:0.6\n"
"      -d, --distance=Num              phasing two variant if distance less than threshold. default:300000\n"
"      -c, --crossBlock=Num            each block tries to connect with next N blocks. default:1\n"
"      -i, --islandBlockLength=Num     phasing across smaller blocks if crossBlock is greater than 1. default:10000\n"
"\n";


static const char* shortopts = "s:b:o:t:r:d:c:i:1:2:3:4:";

enum { OPT_HELP = 1 , DOT_FILE, SV_FILE, IS_ONT, IS_PB};

static const struct option longopts[] = { 
    { "help",                 no_argument,        NULL, OPT_HELP },
    { "dot",                  no_argument,        NULL, DOT_FILE },  
    { "ont",                  no_argument,        NULL, IS_ONT }, 
    { "pb",                   no_argument,        NULL, IS_PB }, 
    { "sv-file",              required_argument,  NULL, SV_FILE },     
    { "reference",            required_argument,  NULL, 'r' },
    { "snp-file",             required_argument,  NULL, 's' },
    { "bam-file",             required_argument,  NULL, 'b' },
    { "out-prefix",           required_argument,  NULL, 'o' },
    { "threads",              required_argument,  NULL, 't' },
    { "distance",             required_argument,  NULL, 'd' },
    { "crossBlock",           required_argument,  NULL, 'c' },
    { "islandBlockLength",    required_argument,  NULL, 'i' },
    { "readsThreshold",       required_argument,  NULL, '1' },
    { "qualityThreshold",     required_argument,  NULL, '2' },
    { "blockReadThreshold",   required_argument,  NULL, '3' },
    { "svReadsThreshold",     required_argument,  NULL, '4' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static int numThreads = 1;
    static int distance = 300000;
    static int crossBlock = 1;
    static int islandBlockLength = 10000;
    static std::string snpFile="";
    static std::string svFile="";
    static std::string bamFile="";
    static std::string fastaFile="";
    static std::string resultPrefix="result";
    static bool generateDot=false;
    static bool isONT=false;
    static bool isPB=false;
    
    double readsThreshold = 0.5;
    double qualityThreshold = 0.7;
    double blockReadThreshold = 1;
    
    double svReadsThreshold = 0.6;
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
        case 'i': arg >> opt::islandBlockLength; break;  

        case '1': arg >> opt::readsThreshold; break; 
        case '2': arg >> opt::qualityThreshold; break; 
        case '3': arg >> opt::blockReadThreshold; break; 
        case '4': arg >> opt::svReadsThreshold; break; 
        
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
            std::cerr<< "File " << opt::snpFile << " not exist.\n\n";
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
            std::cerr<< "File " << opt::bamFile << " not exist.\n\n";
            die = true;
        }
    }
    else{
        std::cerr << SUBPROGRAM ": missing bam file.\n";
        die = true;
    }
    
    if ( opt::numThreads < 1 ){
        std::cerr << SUBPROGRAM " invalid threads. value: " 
                  << opt::numThreads 
                  << "\n please check -t, --threads=Num\n";
        die = true;
    }

    if ( opt::distance < 0 ){
        std::cerr << SUBPROGRAM " invalid distance. value: " 
                  << opt::distance 
                  << "\n please check -d or --distance=Num\n";
        die = true;
    }

    if ( opt::crossBlock < 0 ){
        std::cerr << SUBPROGRAM " invalid crossBlock. value: " 
                  << opt::crossBlock 
                  << "\n please check -c, --crossBlock=Num\n";
        die = true;
    }
    
    if ( opt::islandBlockLength < 0 ){
        std::cerr << SUBPROGRAM " invalid islandBlockLength. value: " 
                  << opt::islandBlockLength 
                  << "\n please check -i, --islandBlockLength=Num\n";
        die = true;
    }

    if ( opt::readsThreshold < 0 || opt::readsThreshold > 1 ){
        std::cerr << SUBPROGRAM " invalid readsThreshold. value: " 
                  << opt::readsThreshold 
                  << "\n please check -1, --readsThreshold=[0~1]\n";
        die = true;
    }
    
    if ( opt::qualityThreshold < 0 || opt::qualityThreshold > 1 ){
        std::cerr << SUBPROGRAM " invalid qualityThreshold. value: " 
                  << opt::qualityThreshold 
                  << "\n please check -2, --qualityThreshold=[0~1]\n";
        die = true;
    }
    
    if ( opt::blockReadThreshold < 0 || opt::blockReadThreshold > 1 ){
        std::cerr << SUBPROGRAM " invalid blockReadThreshold. value: " 
                  << opt::blockReadThreshold 
                  << "\n please check -3, --blockReadThreshold=[0~1]\n";
        die = true;
    }
    
    if ( opt::svReadsThreshold < 0 || opt::svReadsThreshold > 1 ){
        std::cerr << SUBPROGRAM " invalid svReadsThreshold. value: " 
                  << opt::svReadsThreshold 
                  << "\n please check -4, --svReadsThreshold=[0~1]\n";
        die = true;
    }

    if (die)
    {
        std::cerr << "\n" << CORRECT_USAGE_MESSAGE;
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
    ecParams.islandBlockLength=opt::islandBlockLength;
    ecParams.snpFile=opt::snpFile;
    ecParams.svFile=opt::svFile;
    ecParams.bamFile=opt::bamFile;
    ecParams.fastaFile=opt::fastaFile;
    ecParams.resultPrefix=opt::resultPrefix;
    ecParams.generateDot=opt::generateDot;
    ecParams.isONT=opt::isONT;
    ecParams.isPB=opt::isPB;
    
    ecParams.readsThreshold=opt::readsThreshold;
    ecParams.qualityThreshold=opt::qualityThreshold;
    ecParams.blockReadThreshold=opt::blockReadThreshold;
    
    ecParams.svReadsThreshold=opt::svReadsThreshold;
    
    PhasingProcess processor(ecParams);

    return 0;
}