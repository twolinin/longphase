#include "ModCallProcess.h"
#include "ModCallParsingBam.h"

ModCallProcess::ModCallProcess(ModCallParameters params){

    // parsing ref fasta
    std::time_t begin= time(NULL);
    std::cerr<< "reading reference ... ";
    std::map<std::string, std::string> chrString;
    MethFastaParser MethFastaParser(params.fastaFile, chrString);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    bool isFirstchr = true;
    //loop all chromosomes
    for(auto chrIter = chrString.begin(); chrIter != chrString.end(); chrIter++){
        
        // store variant
        std::vector<ReadVariant> fReadVariantVec;
        std::vector<ReadVariant> rReadVariantVec;
        std::vector<ReadVariant> readVariantVec;
        // record hetero metyhl position
        std::map<int,int> *passPosition = new std::map<int,int>;
        
        begin = time(NULL);
        std::string chrName = (*chrIter).first;
        int chrLen = (*chrIter).second.length();
        
        std::cerr<<"parsing contig/chromosome: "<<chrName<<" ... ";
        MethBamParser *methbamparser = new MethBamParser(params, (*chrIter).second);
        
        std::cerr<<"parsing Bam " << " ... ";
        methbamparser->detectMeth(chrName, chrLen, readVariantVec);

        //new new calculate depth
        std::cerr<<"cal depth " << " ... ";
        methbamparser->calculateDepth();
        
        //judge methylation genotype
        std::cerr<<"judge meth " << " ... ";
        methbamparser->judgeMethGenotype(chrName, readVariantVec, fReadVariantVec, rReadVariantVec);
        
        std::cerr<< "run algorithm ... ";
        
        MethylationGraph *fGraph = new MethylationGraph(params);
        MethylationGraph *rGraph = new MethylationGraph(params);
        fGraph->addEdge(fReadVariantVec);
        fGraph->connectResults(chrName, (*passPosition));
        
        rGraph->addEdge(rReadVariantVec);
        rGraph->connectResults(chrName, (*passPosition));
        
        fGraph->destroy();
        rGraph->destroy();
        delete fGraph;
        delete rGraph;

        
        //write to vcf file
        std::cerr<<"write vcf " << " ... ";
        methbamparser->writeResultVCF(chrName, chrString, isFirstchr, (*passPosition));
        isFirstchr = false;
        std::cerr<< difftime(time(NULL), begin) << "s\n";
        
        passPosition->clear();
        
        delete methbamparser;
        delete passPosition;
        
        readVariantVec.clear();
        readVariantVec.shrink_to_fit();
        fReadVariantVec.clear();
        fReadVariantVec.shrink_to_fit();
        rReadVariantVec.clear();
        rReadVariantVec.shrink_to_fit();
    }
}

ModCallProcess::~ModCallProcess(){
};