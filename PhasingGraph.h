#ifndef PHASINGGRAPH_H
#define PHASINGGRAPH_H

#include "Util.h"
#include "ParsingBam.h"
#include "PhasingProcess.h"


typedef std::pair<int, int> PosAllele;
typedef std::map<std::string, int> ReadBase;

class SubEdge{
    
    private:
        int readCount;
        // Edge information. The vector store next pos
        // < next position, read name >
        std::map<int, std::vector<std::string> > refRead;
        std::map<int, std::vector<std::string> > altRead;
        // sum of edge pair quality, pos1 quality + pos2 quality
        // < next position, quality sum >
        std::map<int, int> refQuality;
        std::map<int, int> altQuality;
        // < next position, read count >
        std::map<int, int> refReadCount;
        std::map<int, int> altReadCount;

    public:
        
        SubEdge();
        ~SubEdge();
        
        int getQuality(PosAllele targetPos);
        int getAvgQuality(PosAllele targetPos);
        void addSubEdge(int currentQuality, Variant connectNode, std::string readName);
        std::vector<std::string> showEdge(std::string message);
        std::pair<int,int> BestPair(int targetPos);
        std::vector<std::string> getSupportRead(int breakNodePosition);
};


struct VariantEdge{
    int currPos;
    SubEdge* alt;
    SubEdge* ref;

    VariantEdge(int currPos);
    // node pair 
    std::pair<PosAllele,PosAllele> findBestEdgePair(int targetPos, bool isONT, double diffRatioThreshold, double svReadSimilarRatio, bool conatinSV);
    // number of read of two node. AA and AB combination
    std::pair<int,int> findNumberOfRead(int targetPos);
};


struct BlockRead{
    std::map<std::string,int> readVec;
    
    void recordRead(std::string readName);
};

class VairiantGrpah{
    
    private:
        PhasingParameters *params;
        std::string *ref;
        
        // By default, a Map in C++ is sorted in increasing order based on its key.
        // position, edge
        std::map<int,VariantEdge*> edgeList;
        // record all variant position, include SNP and SV 
        // position, quality
        std::map<int,ReadBase> nodeInfo;
        std::map<int,int> svPosition;
        // phasing result, store final path
        std::map<PosAllele,PosAllele> bestEdgeConnect;
        // the smallest position in the block will be used as the representative of the block
        std::map<int,int>  posAppear;
        std::map<int,int>  blockStart;
        
        // phasing result 
        // PosAllele , block_start    
        std::map<PosAllele,int> bkResult;
        // record each position haplotype
        std::map<PosAllele,int> subNodeHP;
        // store each block and containing positions
        std::map<int,std::vector<int> >  blockVec;

        // homopolymer map
        std::map<int,bool> homopolymerMap;
    
        // disjoint path
        void initialDisjointPath();
        // block phasing
        std::map<std::string,int> getBlockRead(std::pair<int,std::vector<int> > currentBlockVec, std::map<std::string,int> &readQuality, BlockRead &totalRead , int sampleNum);
        bool connectBlockByCommonRead(int nextBlcok, double diffRatioThreshold);
        // connect block by total quality
        void checkTotalQuality();
        
        void checkAverageQuality(double diffRatioThreshold);
        
        void checkDisjointPath(double diffRatioThreshold);
        
        // produce PS tag and determine phased GT tag
        void findResultPath();


    public:
    
        VairiantGrpah(std::string &ref, PhasingParameters &params);
        ~VairiantGrpah();
    
        void addEdge(std::vector<ReadVariant> &readVariant);
        
        void phasingProcess();
        void writingDotFile(std::string dotPrefix);
        void exportResult(std::string chrName, PhasingResult &result);
        int totalNode();
        
};




#endif