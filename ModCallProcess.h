#ifndef MODCALLPROCESS_H
#define MODCALLPROCESS_H
#include "Util.h"


struct ModCallParameters
{
    int numThreads;
    std::string fastaFile;
    std::string resultPrefix;
    std::string methylBamFile;
    float modThreshold;
    float unModThreshold;
    float heterRatio;
    float noiseRatio;
    int connectAdjacent;
    float connectConfidence;
};

class ModCallProcess
{
    public:
        ModCallProcess(ModCallParameters params);
        ~ModCallProcess();
};


#endif