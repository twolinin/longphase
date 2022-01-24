#include "PhasingProcess.h"
#include "PhasingGraph.h"
#include "ParsingBam.h"

PhasingProcess::PhasingProcess(PhasingParameters params)
{
    // load SNP vcf file
    std::time_t begin = time(NULL);
    std::cerr<< "parsing VCF ... ";
    VariantParser snpFile(params.snpFile);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
    
    // load SV vcf file
    begin = time(NULL);
    std::cerr<< "parsing SV VCF ... ";
    SVParser svFile(params.svFile, snpFile);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
 
    // parsing ref fasta 
    begin = time(NULL);
    std::cerr<< "reading reference ... ";
    std::vector<int> last_pos;
    for(auto chr :snpFile.getChrVec()){
        last_pos.push_back(snpFile.getLastSNP(chr));
    }
    FastaParser FastaParser(params.fastaFile , snpFile.getChrVec(), last_pos);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // record all phasing result
    PhasingResult phasingResult;
    // get all detected chromosome
    std::vector<std::string> chrName = snpFile.getChrVec();
 
    // loop all chromosome
    for(std::vector<std::string>::iterator chrIter = chrName.begin(); chrIter != chrName.end() ; chrIter++ ){
        std::cerr<< "parsing contig/chromosome: " << (*chrIter) << " ... ";
        begin = time(NULL);
        
        int lastSNPpos = snpFile.getLastSNP((*chrIter));
        // this chromosome not exist in this file. no variant on this chromosome. 
        if( !snpFile.findChromosome((*chrIter)) || lastSNPpos == -1 ){
            std::cerr<< "skip\n";
            continue;
        }
        
        // store variant
        std::vector<ReadVariant> readVariantVec;

        std::cerr<< "fetch SNP ... ";
        // this method does not store the read information to be used
        BamParser *bamParser = new BamParser((*chrIter), params.bamFile, snpFile, svFile);
        std::string chr_reference = FastaParser.chrString.at(*chrIter);
        bamParser->direct_detect_alleles(lastSNPpos, params, readVariantVec ,chr_reference);
        
        if(params.isONT){
            std::cerr<< "filter SNP ... ";
            snpFile.filterSNP((*chrIter), readVariantVec, chr_reference);
        }

        delete bamParser;
        
        // bam files are partial file or no read support this chromosome's SNP
        if( readVariantVec.size() == 0 ){
            std::cerr<< "skip\n";
            continue;
        }
    
        std::cerr<< "run algorithm ... ";
        
        VairiantGrpah *vGraph = new VairiantGrpah(chr_reference, params);
        // trans read-snp info to edge info
        vGraph->addEdge(readVariantVec);
        // run man algorithm
        vGraph->phasingProcess();
        // push result to phasingResult
        vGraph->exportResult((*chrIter),phasingResult);
        //  generate dot file
        if(params.generateDot)
            vGraph->writingDotFile((*chrIter));
        // free memory
        readVariantVec.clear();
        readVariantVec.shrink_to_fit();
        delete vGraph;
        
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }
    std::cerr<< "writeResult ... ";
    begin = time(NULL);
    snpFile.writeResult(params.resultPrefix, phasingResult);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
    
    if(params.svFile!=""){
        std::cerr<< "write SV Result ... ";
        begin = time(NULL);
        svFile.writeResult(params.resultPrefix, phasingResult);
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }

    return;
};

PhasingProcess::~PhasingProcess(){
    
};

