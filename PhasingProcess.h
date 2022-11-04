#ifndef PHASINGPROCESS_H
#define PHASINGPROCESS_H
#include "Util.h"



struct PhasingParameters
{
    int numThreads;
    int distance;
    int crossSNP;
    //int islandBlockLength;
    std::string snpFile;
    std::string svFile;
    std::string bamFile;
    std::string fastaFile;
    std::string resultPrefix;
    bool generateDot;
    bool isONT;
    bool isPB;
    
    int connectAdjacent;
    int mappingQuality;
    
    double confidentHaplotype;
    double judgeInconsistent;
    int inconsistentThreshold;
    
    double alleleConsistentRatio;
    double maxAlleleRatio;
    
    double readsThreshold;
    //double qualityThreshold;
    //double blockReadThreshold;
    double svReadsThreshold;
    
    std::string version;
    std::string command;
    
    double test1;
    double test2;
};

class PhasingProcess
{

    public:
        PhasingProcess(PhasingParameters params);
        ~PhasingProcess();

};


#endif