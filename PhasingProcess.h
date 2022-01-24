#ifndef PHASINGPROCESS_H
#define PHASINGPROCESS_H
#include "Util.h"



struct PhasingParameters
{
    int numThreads;
    int distance;
    int crossBlock;
    int islandBlockLength;
    std::string snpFile;
    std::string svFile;
    std::string bamFile;
    std::string fastaFile;
    std::string resultPrefix;
    bool generateDot;
    bool isONT;
    bool isPB;
    
    double readsThreshold;
    double qualityThreshold;
    double blockReadThreshold;
    double svReadsThreshold;
};

class PhasingProcess
{

    public:
        PhasingProcess(PhasingParameters params);
        ~PhasingProcess();

};


#endif