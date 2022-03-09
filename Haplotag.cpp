#include "Haplotag.h"
#include "HaplotagProcess.h"
#include "Util.h"
#include <getopt.h>


#define SUBPROGRAM "haplotag"

static const char *CORRECT_USAGE_MESSAGE =
"Usage: "  " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"      --help                          display this help and exit.\n\n"
"require arguments:\n"
"      -s, --snp-file=NAME             input SNP vcf file.\n"
"      -b, --bam-file=NAME             input bam file.\n"
"optional arguments:\n"
"      --tagSupplementary              tag supplementary alignment. default:false\n"
"      -q, --qualityThreshold=Num      not tag alignment if the mapping quality less than threshold. default:0\n"
"      -p, --percentageThreshold=Num   the alignment will be tagged according to the haplotype corresponding to most alleles.\n"
"                                      if the alignment has no obvious corresponding haplotype, it will not be tagged. default:0.6\n"
"      -t, --threads=Num               number of thread. default:1\n"
"      -o, --out-prefix=NAME           prefix of phasing result. default:result\n";



static const char* shortopts = "s:b:o:t:q:p:";

enum { OPT_HELP = 1, TAG_SUP};

static const struct option longopts[] = { 
    { "help",                 no_argument,        NULL, OPT_HELP },
    { "snp-file",             required_argument,  NULL, 's' },
    { "bam-file",             required_argument,  NULL, 'b' },
    { "tagSupplementary",     no_argument,        NULL, TAG_SUP },
    { "out-prefix",           required_argument,  NULL, 'o' },
    { "threads",              required_argument,  NULL, 't' },
    { "qualityThreshold",     required_argument,  NULL, 'q' },
    { "percentageThreshold",  required_argument,  NULL, 'p' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static int numThreads = 1;
    static int qualityThreshold = 20;
    static double percentageThreshold = 0.6;
    static std::string snpFile="";
    static std::string bamFile="";
    static std::string resultPrefix="result";
    static bool tagSupplementary = false;
}

void HaplotagOptions(int argc, char** argv)
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
        case 'q': arg >> opt::qualityThreshold; break;
        case 'p': arg >> opt::percentageThreshold; break;
        case TAG_SUP:
             opt::tagSupplementary = true;
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
                  << "\nplease check -t, --threads=Num\n";
        die = true;
    }
    
    if ( opt::percentageThreshold > 1 || opt::percentageThreshold < 0 ){
        std::cerr << SUBPROGRAM " invalid percentage threshold. value: " 
                  << opt::percentageThreshold
                  << "\nthis value need: 0~1, please check -p, --percentageThreshold=Num\n";
        die = true;
    }
    
    if (die)
    {
        std::cerr << "\n" << CORRECT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

}

int HaplotagMain(int argc, char** argv)
{
    HaplotagParameters ecParams;
    // set parameters
    HaplotagOptions(argc, argv);

    ecParams.numThreads=opt::numThreads;
    ecParams.qualityThreshold=opt::qualityThreshold;
    ecParams.snpFile=opt::snpFile;
    ecParams.bamFile=opt::bamFile;
    ecParams.resultPrefix=opt::resultPrefix;
    ecParams.tagSupplementary=opt::tagSupplementary;
    ecParams.percentageThreshold=opt::percentageThreshold;

    HaplotagProcess processor(ecParams);

    return 0;
}