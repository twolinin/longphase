#ifndef MODCALLPARSINGBAM_H
#define MODCALLPARSINGBAM_H

#include "Util.h"
#include "ModCallProcess.h"
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/kbitset.h>
#include <htslib/thread_pool.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kstring.h>

#include "ParsingBam.h"
#include "PhasingGraph.h"
#include "PhasingProcess.h"


struct MethPosInfo{
    MethPosInfo():methreadcnt(0), noisereadcnt(0),canonreadcnt(0), depth(0),heterstatus(""),strand(-1){}

    int methreadcnt;
    int noisereadcnt;
    int canonreadcnt;
    int depth;
    std::string heterstatus;
    int strand; //0: + , 1: -
    std::vector<std::string> modReadVec;
    std::vector<std::string> nonModReadVec;
};

// The ReferenceChromosome struct is used to store information about a reference fasta.
struct ReferenceChromosome {
    std::string name;      
    std::string sequence;  
    int length;            

    ReferenceChromosome(const std::string& name, const std::string& sequence, int length)
        : name(name), sequence(sequence), length(length) {}
};

//FastaParser: to get each chromosome name and lastpos(for chunksize)
class MethFastaParser{
    private:
        // file name
        std::string fastaFile ;

    public:
        MethFastaParser(std::string fastaFile, std::vector<ReferenceChromosome> &chrInfo);
        ~MethFastaParser();
};

struct MethPosProb{
    MethPosProb(int position, int prob):
    position(position),
    prob(prob){};
    
    int position;
    int prob;
};

struct AlignmentMethExtend: Alignment{    
    //hts_base_mod_state *m;
    //store the read position and probability
    std::vector<MethPosProb> queryMethVec;
};

class MethBamParser{
    private:
        ModCallParameters *params;
        std::string *refString;
        std::string chrName;
        int refstartpos;
        
        std::map<int , MethPosInfo> *chrMethMap;
        //<pos, <positive strand cnt, negative strand cnt>>
        std::map<int, std::pair<int,int>> *readStartEndMap; 
        
        // get methylation tag from BAM file
        std::string get_aux_tag(const bam1_t *aln, const char tag[2]);
        void parse_CIGAR(const  bam_hdr_t &bamHdr,const bam1_t &aln, std::vector<ReadVariant> &readVariantVec);
        
        // record methylation position and probability from methylation tag
        void getmeth(AlignmentMethExtend &align);
        
    public:
        MethBamParser(ModCallParameters &params, std::string &refString);
        ~MethBamParser();
        void detectMeth(std::string chrName, int chr_len, int numThreads, std::vector<ReadVariant> &readVariantVec);
        void exportResult(std::string chrName, std::string chrSquence, int chrLen , std::map<int,int> &passPosition, std::ostringstream &methResult);
        void judgeMethGenotype(std::string chrName, std::vector<ReadVariant> &readVariantVec, std::vector<ReadVariant> &fReadVariantVec, std::vector<ReadVariant> &rReadVariantVec );
        void calculateDepth();
        
};

void writeResultVCF(ModCallParameters &params, std::vector<ReferenceChromosome> &chrInfo, std::map<std::string,std::ostringstream> &chrResult);

class MethylationGraph{
    private:
        ModCallParameters *params;
        std::vector<ReadVariant> *readVariant;
        
        // modifications
        // position, quality
        std::map<int,ReadBaseMap*> *forwardModNode;
        std::map<int,ReadBaseMap*> *reverseModNode;
        
        // By default, a Map in C++ is sorted in increasing order based on its key.
        // position, edge
        std::map<int,VariantEdge*> *edgeList;
        // record all variant position, include SNP and SV 
        // position, quality
        std::map<int,ReadBaseMap*> *nodeInfo;
        
    public:
    
        MethylationGraph(ModCallParameters &params);
        ~MethylationGraph();
        
        void addEdge(std::vector<ReadVariant> &readVariant);
        void connectResults(std::string chrName, std::map<int,int> &passPosition);
        
        void destroy();
};

#endif