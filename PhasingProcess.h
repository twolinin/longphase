#ifndef PHASINGPROCESS_H
#define PHASINGPROCESS_H
#include "Util.h"



struct PhasingParameters
{
    int numThreads;
    int distance;
    std::string snpFile;
    std::string svFile;
    std::vector<std::string> bamFile;
    std::string modFile="";
    std::string fastaFile;
    std::string resultPrefix;
    bool generateDot;
    bool isONT;
    bool isPB;
    bool phaseIndel;
    
    int connectAdjacent;
    int mappingQuality;
    double mismatchRate;
    
    int baseQuality;
    double edgeWeight;
    //double falseWeight;
    //int coverageThreshold;
    
    double snpConfidence;
    double readConfidence;
    
    double edgeThreshold;
    
    std::string version;
    std::string command;
};

class PhasingProcess
{

    public:
        PhasingProcess(PhasingParameters params);
        ~PhasingProcess();

};


#endif
