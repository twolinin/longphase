#ifndef PHASINGGRAPH_H
#define PHASINGGRAPH_H

#include "Util.h"
#include "ParsingBam.h"
#include "PhasingProcess.h"


typedef std::pair<int, int> PosAllele;
typedef std::map<std::string, int> ReadBaseMap;
using PosVec = std::vector<int>;

class Clip{
    private:
        std::string chr;
        PosVec CNVtoLOHInterval(SnpParser &snpMap);
        //std::map<int,int> cnvMap;
        struct cnvState{
            bool push;
            bool slowUp;
            bool slowDown;
            int currCount;
            int rejectCount;
            int pullDownCount;
            int slowDownCount;
            int candidateStartPos;
            int candidateEndPos;
            
            void reset(){
                push = false;
                slowUp = false;
                slowDown = false;
                currCount = 0;
                rejectCount = 0;
                pullDownCount = 0;
                slowDownCount = 0;
                candidateStartPos = -1;
                candidateEndPos = -1;
            }
        };
        cnvState state;
        void updateThreshold(int upCount);

    public:
        Clip(std::string &chr, ClipCount &inClipCount);
        ~Clip();
        //const std::map<int, int>& getCnvMap() const{return cnvMap;}
        std::map<int,int> cnvMap;
        std::vector<std::pair<int, int>> cnvVec;
        void getCNVInterval(ClipCount &clipCount);
        PosVec detectLOH(SnpParser &snpMap);
        std::string getChr() const { return chr;}
};

class SubEdge{
    
    private:
        int readCount;
        // Edge information. The vector store next pos
        // < next position, read name >
        std::map<int, std::vector<std::string> > *refRead;
        std::map<int, std::vector<std::string> > *altRead;
        // sum of edge pair quality, pos1 quality + pos2 quality
        // < next position, quality sum >
        std::map<int, int> *refQuality;
        std::map<int, int> *altQuality;
        // < next position, read count >
        std::map<int, float> *refReadCount;
        std::map<int, float> *altReadCount;

    public:
        
        SubEdge();
        ~SubEdge();
        
        void destroy();
        
        void addSubEdge(int currentQuality, Variant connectNode, std::string readName, int baseQuality, double edgeWeight);
        std::pair<float,float> BestPair(int targetPos);
        float getRefReadCount(int targetPos);
        float getAltReadCount(int targetPos);        
        
        std::vector<std::string> showEdge(std::string message);
        std::vector<std::pair<int,int>> getConnectPos();
 
        int getQuality(PosAllele targetPos);
        int getAvgQuality(PosAllele targetPos);

};

//use to store the voting info from the previous variants
struct VoteResult{
    int Pos;		//who votes
    float para;		//rr+aa
    float cross;	//ra+ar
    float weight;	//how much weight
    int hap;		//which haplotype 
    double ESR;		//similarity of para and cross

    VoteResult( int currPos, float weight ) ;
};

struct VariantEdge{
    std::string chrName;
    int currPos;
    SubEdge* alt;
    SubEdge* ref;
    VariantEdge(int currPos);
    // node pair 
    std::pair<PosAllele,PosAllele> findBestEdgePair(int targetPos, bool isONT, double diffRatioThreshold, bool debug, std::map<int,int> &variantType, VoteResult &vote);
    // number of read of two node. AA and AB combination
    std::pair<float,float> findNumberOfRead(int targetPos);
};


struct BlockRead{
    std::map<std::string,int> readVec;
    
    void recordRead(std::string readName);
};

struct EdgeResult{
    int rr;
    int ra;
    int ar;
    int aa;
};

struct CnvStatistics{
    std::map<int, double> missRateMap;
    std::map<int, std::map<int, std::vector<int>>> cnvReadMmrate;
    std::map<int, std::map<int, double>> AvgCnvReadMiss;
};

class VairiantGraph{
    
    private:
        PhasingParameters *params;
        std::string *ref;
        std::vector<std::string> dotResult;
        std::vector<ReadVariant> *readVariant;
        
        // By default, a Map in C++ is sorted in increasing order based on its key.
        // position, edge
        std::map<int,VariantEdge*> *edgeList;

        // Each position will record the included reads and their corresponding base qualities.
        // position, < read name, quality>
        std::map<int,ReadBaseMap*> *totalVariantInfo;
        // position, type < 0=SNP 1=SV 2=MOD 3=INDEL >
        std::map<int,int> *variantType;

        std::pair<float,float> Onelongcase( std::vector<VoteResult> vote ) ;

        // phasing result     
        // PosAllele , block_start    
        std::map<PosAllele,int> *bkResult;
        // record each position haplotype
        std::map<PosAllele,int> *subNodeHP;
        // store phased read and read's haplotype
        std::map<std::string,int> *readHpMap;

        // produce PS tag and determine phased GT tag
        void storeResultPath();
        
        void readCorrection();

        void edgeConnectResult();

        void calculateCnvMismatchRate(std::vector<ReadVariant>& in_readVariant, Clip &clip);

        void aggregateCnvReadMismatchRate(const std::vector<ReadVariant>& in_readVariant, const Clip &clip, std::map<int, std::map<int, std::vector<int>>>& cnvReadMmrate);

        void calculateAverageMismatchRate(const Clip& clip, const std::map<int, std::map<int, std::vector<int>>>& cnvReadMmrate, std::map<int, double>& missRateMap);

        void filterHighMismatchVariants(std::vector<ReadVariant>& in_readVariant, const Clip& clip, const std::map<int, double>& missRateMap);

        bool isPositionInRange(int position, int start, int end);
        

    public:
    
        VairiantGraph(std::string &ref, PhasingParameters &params);
        ~VairiantGraph();
    
        void addEdge(std::vector<ReadVariant> &in_readVariant, Clip &clip);
        
        void phasingProcess();
        void writingDotFile(std::string dotPrefix);
        std::map<std::string,int>* getReadHP();
        void exportResult(std::string chrName, PhasingResult &result);
        int totalNode();

        void destroy();
        
};



#endif
