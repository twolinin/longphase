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
"      -r, --reference=NAME            reference fasta.\n"
"optional arguments:\n"
"      --tagSupplementary              tag supplementary alignment. default:false\n"
"      --sv-file=NAME                  input phased SV vcf file.\n"
"      --mod-file=NAME                 input a modified VCF file (produced by longphase modcall and processed by longphase phase).\n"
"      -q, --qualityThreshold=Num      not tag alignment if the mapping quality less than threshold. default:1\n"
"      -p, --percentageThreshold=Num   the alignment will be tagged according to the haplotype corresponding to most alleles.\n"
"                                      if the alignment has no obvious corresponding haplotype, it will not be tagged. default:0.6\n"
"      -t, --threads=Num               number of thread. default:1\n"
"      -o, --out-prefix=NAME           prefix of phasing result. default:result\n"
"      --region=REGION                 tagging include only reads/variants overlapping those regions. default:""(all regions)"
"      --log                           an additional log file records the result of each read. default:false\n";


static const char* shortopts = "s:b:o:t:q:p:r:";

enum { OPT_HELP = 1, TAG_SUP, SV_FILE, REGION, LOG, MOD_FILE};

static const struct option longopts[] = { 
    { "help",                 no_argument,        NULL, OPT_HELP },
    { "snp-file",             required_argument,  NULL, 's' },
    { "bam-file",             required_argument,  NULL, 'b' },
    { "reference",            required_argument,  NULL, 'r' },
    { "tagSupplementary",     no_argument,        NULL, TAG_SUP },
    { "sv-file",              required_argument,  NULL, SV_FILE },
    { "mod-file",             required_argument,  NULL, MOD_FILE },
    { "out-prefix",           required_argument,  NULL, 'o' },
    { "threads",              required_argument,  NULL, 't' },
    { "qualityThreshold",     required_argument,  NULL, 'q' },
    { "percentageThreshold",  required_argument,  NULL, 'p' },
    { "region",               required_argument,  NULL, REGION },
    { "log",                  no_argument,        NULL, LOG },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static int numThreads = 1;
    static int qualityThreshold = 1;
    static double percentageThreshold = 0.6;
    static std::string snpFile="";
    static std::string svFile="";
    static std::string modFile="";
    static std::string bamFile="";
    static std::string fastaFile="";
    static std::string resultPrefix="result";
    static std::string region="";
    static bool tagSupplementary = false;
    static bool writeReadLog = false;
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
        case 'r': arg >> opt::fastaFile; break; 
        case 'o': arg >> opt::resultPrefix; break;
        case 'q': arg >> opt::qualityThreshold; break;
        case 'p': arg >> opt::percentageThreshold; break;
        case SV_FILE: arg >> opt::svFile; break;
        case MOD_FILE: arg >> opt::modFile; break;
        case TAG_SUP: opt::tagSupplementary = true; break;
        case REGION: arg >> opt::region; break;
        case LOG: opt::writeReadLog = true; break;
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
    
    if( opt::svFile != "")
    {
        std::ifstream openFile( opt::svFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cerr<< "File " << opt::svFile << " not exist.\n\n";
            die = true;
        }
    }
    
    if( opt::modFile != "")
    {
        std::ifstream openFile( opt::modFile.c_str() );
        if( !openFile.is_open() )
        {
            std::cerr<< "File " << opt::modFile << " not exist.\n\n";
            die = true;
        }
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
    
    if( opt::fastaFile != "")
    {
        std::ifstream openFile( opt::snpFile.c_str() );
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
    ecParams.svFile=opt::svFile;
    ecParams.modFile=opt::modFile;
    ecParams.bamFile=opt::bamFile;
    ecParams.fastaFile=opt::fastaFile;
    ecParams.resultPrefix=opt::resultPrefix;
    ecParams.tagSupplementary=opt::tagSupplementary;
    ecParams.percentageThreshold=opt::percentageThreshold;
    ecParams.region=opt::region;
    ecParams.writeReadLog=opt::writeReadLog;

    HaplotagProcess processor(ecParams);

    return 0;
}
