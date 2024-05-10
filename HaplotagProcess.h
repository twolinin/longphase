#ifndef HAPLOTAGPROCESS_H
#define HAPLOTAGPROCESS_H

#include "Util.h"
#include "ParsingBam.h"
#include <cmath>

struct HaplotagParameters
{
    int numThreads;
    int qualityThreshold;
    
    double percentageThreshold;
    
    std::string snpFile;
    std::string svFile;
    std::string modFile;
    std::string bamFile;
    std::string fastaFile;
    std::string resultPrefix;
    std::string region;
    std::string command;
    std::string version;
    
    bool tagSupplementary;
    bool writeReadLog;
};

class HaplotagProcess
{
    void variantParser(std::string variantFile);
    void compressParser(std::string &variantFile);
    void unCompressParser(std::string &variantFile);
    void parserProcess(std::string &input);
    
    void tagRead(HaplotagParameters &params);
    
    std::vector<std::string> chrVec;
    std::map<std::string, int> chrLength;
    
    // chr, variant position (0-base), allele haplotype
    std::map<std::string, std::map<int, RefAlt> > chrVariant;
    // chr, variant position (0-base), phased set
    std::map<std::string, int > psIndex;
    std::map<std::string, std::map<int, int> > chrVariantPS;
    
    std::map<std::string, std::map<int, std::string> > chrVariantHP1;
    std::map<std::string, std::map<int, std::string> > chrVariantHP2;

    std::map<int, RefAlt> currentChrVariants;
    std::map<int, RefAlt>::iterator firstVariantIter;
    // The number of SVs occurring on different haplotypes in a read
    std::map<std::string, std::map<int, int> > readSVHapCount;

    void initFlag(bam1_t *aln, std::string flag);
    
    int judgeHaplotype(const  bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, const std::string &ref_string);
    
    int totalAlignment;
    int totalSupplementary;
    int totalSecondary;
    int totalUnmapped;
    int totalTagCount;
    int totalUnTagCount;
    
    std::time_t processBegin;
    bool integerPS;
    bool parseSnpFile;
    bool parseSVFile;
    bool parseMODFile;
    
    public:
        HaplotagProcess(HaplotagParameters params);
        ~HaplotagProcess();

};


#endif