#ifndef PHASINGGRAPH_H
#define PHASINGGRAPH_H

#include "Util.h"
#include "ParsingBam.h"
#include "PhasingProcess.h"
#include <set>
#include <array>
#include <algorithm>
#include <cstdint>


typedef std::pair<int, int> PosAllele;
typedef std::map<std::string, int> ReadBaseMap;
// [haplotype 0/1][position][allele 0/1] -> read count; haplotype and allele index fixed arrays
// [haplotype][variant index] -> < ref count, alt count >
using HpAlleleCountMap = std::array<std::vector<std::array<double, 2>>, 2>;
using PhasedBlocks = std::map<int, std::vector<int>>;

class Clip{
    private:
        std::string chr;

        static const int CNV_AREA_SIZE             = 30000;
        static const int CNV_PUSH_THRESHOLD        = 5;
        static const int CNV_EXTEND_THRESHOLD      = 30;
        static const int CNV_MAX_LENGTH            = 200000;
        static const int CNV_SLOWUP_COMMIT_THRESHOLD = 20;

        struct cnvState{
            cnvState():push(false),slowUp(false),slowDown(false),currCount(0),rejectCount(0),pullDownCount(0),slowDownCount(0),candidateStartPos(-1),candidateEndPos(-1){}
            
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
        void handleIdleState(int upCount, int downCount, int pos);
        void handleActiveState(int upCount, int downCount, int pos);
        void handleSlowUpState(int upCount, int downCount, int pos);

    public:
        Clip(std::string &chr, ClipCount &inClipCount);
        ~Clip();
        std::vector<std::pair<int, int>> cnvVec;
        void getCNVInterval(ClipCount &clipCount);
};

class SubEdge{
    
    private:
        int readCount;
        // < next position, read count >, sorted by position. A sorted array instead of std::map:
        // each table holds a few positions that mostly arrive in increasing order, so lookups
        // stay in one or two cache lines instead of following tree pointers.
        typedef std::vector<std::pair<int, float>> CountTable;
        CountTable refReadCount;
        CountTable altReadCount;
        // index of the last position touched in each table; the next lookup starts there
        size_t refHint;
        size_t altHint;

        static float &findOrInsert(CountTable &table, size_t &hint, int pos);
        static float findCount(const CountTable &table, int pos);

    public:
        
        SubEdge();
        ~SubEdge();
        
        void destroy();
        
        void addSubEdge(int currentQuality, const Variant &connectNode, const std::string &readName, int baseQuality, double edgeWeight);
        std::pair<float,float> BestPair(int targetPos);
        float getRefReadCount(int targetPos);
        float getAltReadCount(int targetPos);        
        
        std::vector<std::string> showEdge(std::string message);
        std::vector<std::pair<int,int>> getConnectPos();

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
        std::string *chrName;
        PhasingParameters *params;
        std::string *ref;
        std::string dotResult;
        // graph edges for the GNN
        std::vector<DotEdge> gnnEdges;
        std::vector<ReadVariant> *readVariant;
        
        // Variant positions in genomic order; every per-variant array below is indexed
        // like this vector.
        std::vector<int> variantPositions;
        // read support from each variant to the next connectAdjacent variants (the scan never
        // looks further): edgeCounts[index * connectAdjacent + distance - 1] holds
        // < ref->ref, ref->alt, alt->ref, alt->alt > weights
        std::vector<std::array<float, 4> > edgeCounts;
        // the variant starts an edge (some read has a variant after it)
        std::vector<char> hasEdge;
        // type < 0=SNP 1=SV 2=MOD 3=INDEL 4=danger INDEL >
        std::vector<int8_t> variantType;

        std::pair<float,float> Onelongcase( const std::vector<VoteResult> &vote ) ;

        // node pair between currPos and targetPos from the four allele-to-allele read weights
        std::pair<PosAllele,PosAllele> findBestEdgePair(int currPos, int targetPos, float rr, float ra, float ar, float aa, double edgeThreshold, bool debug, int currType, int targetType, VoteResult &vote);

        // phasing result, per variant index * 2 + (0 ref, 1 alt):
        // block start + 1 (0 = not in a block)
        std::vector<int> bkResult;
        // haplotype of each allele
        std::vector<int8_t> subNodeHP;
        // per-variant vote diagnostics recorded during
        // scanVariantsAndBuildBlocks() and consumed in exportResult():
        // weighted HP1 / HP2 vote counts and the Shannon entropy of (h1, h2) in bits (0..1)
        std::vector<float> h1weight;
        std::vector<float> h2weight;
        std::vector<float> variantEntropy;
        std::vector<char> hasEntropy;

        // index of a variant position, -1 if absent
        int variantIndex(int position) const;

        // produce PS tag and determine phased GT tag
        void storeResultPath();
        
        void readCorrection();

        void accumulateReadHaplotypeEvidence(const Variant& variant, double& refCount, double& altCount);

        void assignReadHaplotypesAndCollectAlleleCounts(HpAlleleCountMap& hpAlleleCountMap);

        void reassignVariantHaplotypes(const HpAlleleCountMap& hpAlleleCountMap);

        void edgeConnectResult();

        void scanVariantsAndBuildBlocks(
            std::vector<int>& hpResult,
            PhasedBlocks& phasedBlocks,
            std::string& outDotResult);

        void materializeBlockResults(
            const std::vector<int>& hpResult,
            const PhasedBlocks& phasedBlocks);

        void calculateCnvMismatchRate(std::vector<ReadVariant>& in_readVariant, Clip &clip);

        void aggregateCnvReadMismatchRate(const std::vector<ReadVariant>& in_readVariant, const Clip &clip, std::map<int, std::map<int, std::vector<int>>>& cnvReadMmrate);

        void calculateAverageMismatchRate(const Clip& clip, const std::map<int, std::map<int, std::vector<int>>>& cnvReadMmrate, std::map<int, double>& missRateMap);

        void filterHighMismatchVariants(std::vector<ReadVariant>& in_readVariant, const Clip& clip, const std::map<int, double>& missRateMap);

        void filterOverlappingAlignments(std::vector<ReadVariant>& in_readVariant);

        void applyCnvFilter(std::vector<ReadVariant>& in_readVariant, Clip& clip);

        void buildVariantGraph(std::vector<ReadVariant>& in_readVariant);

        bool isPositionInRange(int position, int start, int end);
        

    public:
    
        VairiantGraph(std::string &ref, PhasingParameters &params, std::string &chrName);
        ~VairiantGraph();
    
        void addEdge(std::vector<ReadVariant> &in_readVariant, Clip &clip);
        
        void phasingProcess();
        void writingDotFile(std::string dotPrefix);
        void exportGnnEdges(std::vector<DotEdge> &out);
        void exportResult(std::string chrName, PhasingResult &result);
        int totalNode();

        void destroy();
        
};



#endif
