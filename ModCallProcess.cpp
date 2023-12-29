#include "ModCallProcess.h"
#include "ModCallParsingBam.h"

ModCallProcess::ModCallProcess(ModCallParameters params){

    // parsing ref fasta
    std::time_t begin= time(NULL);
    std::cerr<< "reading reference ... ";
    std::vector<ReferenceChromosome> chrInfo;
    MethFastaParser MethFastaParser(params.fastaFile, chrInfo);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
    
    // record all ModCall result
    std::map<std::string,std::ostringstream> chrModCallResult;
    // Initialize an empty map in chrModCallResult to store all ModCall results.
    // This is done to prevent issues with multi-threading, by defining an empty map first.
    for (auto chrIter = chrInfo.begin(); chrIter != chrInfo.end(); chrIter++) {
        chrModCallResult[chrIter->name] = std::ostringstream();
    }

    // set chrNumThreads and bamParserNumThreads based on parameters
    int chrNumThreads,bamParserNumThreads;
    setModcallNumThreads(params.numThreads, chrNumThreads, bamParserNumThreads);
    // setNumThreads(chrInfo.size(), params.numThreads, chrNumThreads, bamParserNumThreads);
    begin = time(NULL);

    //loop all chromosomes
    #pragma omp parallel for schedule(dynamic) num_threads(chrNumThreads)
    for(auto chrIter = chrInfo.begin(); chrIter != chrInfo.end(); ++chrIter) {
        std::time_t chrbegin = time(NULL);

        // store variant
        std::vector<ReadVariant> fReadVariantVec;
        std::vector<ReadVariant> rReadVariantVec;
        std::vector<ReadVariant> readVariantVec;
        // record hetero metyhl position
        std::map<int,int> *passPosition = new std::map<int,int>;

        std::string chrName = chrIter->name;
        std::string chrSeq = chrIter->sequence;
        int chrLen = chrIter->length;
        
        MethBamParser *methbamparser = new MethBamParser(params, chrSeq);
        
        methbamparser->detectMeth(chrName, chrLen, bamParserNumThreads, readVariantVec);

        //new new calculate depth
        methbamparser->calculateDepth();
        
        //judge methylation genotype
        methbamparser->judgeMethGenotype(chrName, readVariantVec, fReadVariantVec, rReadVariantVec);

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
        
        //push result to ModCallResult
        methbamparser->exportResult(chrName, chrSeq, chrLen, (*passPosition), chrModCallResult[chrName]);
        
        passPosition->clear();
        
        delete methbamparser;
        delete passPosition;
        
        readVariantVec.clear();
        readVariantVec.shrink_to_fit();
        fReadVariantVec.clear();
        fReadVariantVec.shrink_to_fit();
        rReadVariantVec.clear();
        rReadVariantVec.shrink_to_fit();

        std::cerr<< "(" << chrName << "," << difftime(time(NULL), chrbegin) << "s)";
    }
    std::cerr<< "\nmodcall total:  " << difftime(time(NULL), begin) << "s\n";

    begin= time(NULL);
    std::cerr<<"write vcf " << " ... ";
    //write to vcf file
    writeResultVCF(params, chrInfo, chrModCallResult);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
}

ModCallProcess::~ModCallProcess(){
};