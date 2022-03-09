#ifndef HAPLOTAGPROCESS_H
#define HAPLOTAGPROCESS_H

#include "Util.h"
#include "ParsingBam.h"


struct HaplotagParameters
{
    int numThreads;
    int qualityThreshold;
    
    double percentageThreshold;
    
    std::string snpFile;
    std::string bamFile;
    std::string resultPrefix;
    
    bool tagSupplementary;

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

    std::map<int, RefAlt> currentVariants;
    std::map<int, RefAlt>::iterator firstVariantIter;

    void initFlag(bam1_t *aln, std::string flag);
    
    int judgeHaplotype(Alignment align, std::string chrName, double percentageThreshold);
    
    int totalAlignment;
    int totalSupplementary;
    int totalSecondary;
    int totalUnmapped;
    int totalTagCuonnt;
    int totalUnTagCuonnt;
    
    std::time_t processBegin;
    bool integerPS;
    
    public:
        HaplotagProcess(HaplotagParameters params);
        ~HaplotagProcess();

};


#endif