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
        int getRefReadCount(int targetPos);
        int getAltReadCount(int targetPos);
        
        void addSubEdge(int currentQuality, Variant connectNode, std::string readName);
        std::vector<std::string> showEdge(std::string message);
        std::pair<int,int> BestPair(int targetPos);
        void getReadVariant(std::map< std::string,std::map<int,int> > &readVariantMap);
};


struct VariantEdge{
    int currPos;
    SubEdge* alt;
    SubEdge* ref;

    VariantEdge(int currPos);
    // node pair 
    std::pair<PosAllele,PosAllele> findBestEdgePair(int targetPos, bool isONT, double diffRatioThreshold, bool repair);
    // number of read of two node. AA and AB combination
    std::pair<int,int> findNumberOfRead(int targetPos);
};


struct BlockRead{
    std::map<std::string,int> readVec;
    
    void recordRead(std::string readName);
};

class VairiantGraph{
    
    private:
        PhasingParameters *params;
        std::string *ref;
        std::vector<std::string> dotResult;
        std::vector<ReadVariant> *readVariant;
        
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
        // current snp, result haplotype (1 or 2)
        std::map<int,int> hpResult;        
        // PosAllele , block_start    
        std::map<PosAllele,int> bkResult;
        // record each position haplotype
        std::map<PosAllele,int> subNodeHP;
        // store each block and containing positions
        std::map<int,std::vector<int> >  blockVec;
        
        // block phasing
        std::vector<ReadVariant> getBlockRead(std::pair<int,std::vector<int> > currentBlockVec, BlockRead &totalRead , int sampleNum);
        bool blockPhaseing();

        // produce PS tag and determine phased GT tag
        void storeResultPath();
        
        void readCorrection();

        void edgeConnectResult();

    public:
    
        VairiantGraph(std::string &ref, PhasingParameters &params);
        ~VairiantGraph();
    
        void addEdge(std::vector<ReadVariant> &readVariant);
        
        void phasingProcess();
        void writingDotFile(std::string dotPrefix);
        void exportResult(std::string chrName, PhasingResult &result);
        int totalNode();
        
};




#endif