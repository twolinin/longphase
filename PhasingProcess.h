#ifndef PHASINGPROCESS_H
#define PHASINGPROCESS_H
#include "Util.h"



struct PhasingParameters
{
    int numThreads;
    int distance;
    int crossBlock;
    std::string snpFile;
    std::string svFile;
    std::string bamFile;
    std::string fastaFile;
    std::string resultPrefix;
    bool generateDot;
    bool isONT;
    bool isPB;
    
    double readSimilarRatio;
    double quaSimilarRatio;
    double blockSimilarRatio;
};

class PhasingProcess
{

    public:
        PhasingProcess(PhasingParameters params);
        ~PhasingProcess();

};


#endif