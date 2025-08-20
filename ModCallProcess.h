#ifndef MODCALLPROCESS_H
#define MODCALLPROCESS_H
#include "Util.h"


struct ModCallParameters
{
    int numThreads;
    std::string fastaFile;
    std::string snpFile;
    std::string resultPrefix;
    std::vector<std::string> bamFileVec;
    float modThreshold;
    float unModThreshold;
    float heterRatio;
    float noiseRatio;
    int connectAdjacent;
    float connectConfidence;
    int iterCount;
    
    std::string version;
    std::string command;
    
    bool outputAllMod;
};

class ModCallProcess
{
    public:
        ModCallProcess(ModCallParameters params);
        ~ModCallProcess();
};


#endif