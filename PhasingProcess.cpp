#include "PhasingProcess.h"
#include "PhasingGraph.h"
#include "ParsingBam.h"

PhasingProcess::PhasingProcess(PhasingParameters params)
{
    std::cerr<< "LongPhase Ver " << params.version << "\n";
    std::cerr<< "\n";
    std::cerr<< "--- File Parameter --- \n";
    std::cerr<< "SNP File           : " << params.snpFile      << "\n";
    std::cerr<< "SV  File           : " << params.svFile       << "\n";
    std::cerr<< "MOD File           : " << params.modFile      << "\n";
    std::cerr<< "REF File           : " << params.fastaFile    << "\n";
    std::cerr<< "Output Prefix      : " << params.resultPrefix << "\n";
    std::cerr<< "Number of Threads  : " << params.numThreads          << "\n";
    std::cerr<< "Generate Dot       : " << ( params.generateDot ? "True" : "False" ) << "\n";
    std::cerr<< "BAM File           : ";
    for( auto file : params.bamFile){
        std::cerr<< file <<" " ;   
    }
    std::cerr << "\n";
    
    std::cerr<< "\n";
    std::cerr<< "--- Phasing Parameter --- \n";
    std::cerr<< "Seq Platform       : " << ( params.isONT ? "ONT" : "PB" ) << "\n";
    std::cerr<< "Phase Indel        : " << ( params.phaseIndel ? "True" : "False" )  << "\n";
    std::cerr<< "Distance Threshold : " << params.distance        << "\n";
    std::cerr<< "Connect Adjacent   : " << params.connectAdjacent << "\n";
    std::cerr<< "Edge Threshold     : " << params.edgeThreshold   << "\n";
    std::cerr<< "Overlap Threshold  : " << params.overlapThreshold   << "\n";
    std::cerr<< "Mapping Quality    : " << params.mappingQuality  << "\n";
    std::cerr<< "Mismatch Rate      : " << params.mismatchRate  << "\n";
    std::cerr<< "Variant Confidence : " << params.snpConfidence   << "\n";
    std::cerr<< "ReadTag Confidence : " << params.readConfidence  << "\n";
    if (!params.svFile.empty()) {
        std::cerr<< "SV Windowsize      : " << params.svWindow << "\n";
        std::cerr<< "SV Threshold       : " << params.svThreshold << "\n";
    }
    std::cerr<< "\n";
    
    std::time_t processBegin = time(NULL);
        
    // load SNP vcf file
    std::time_t begin = time(NULL);
    std::cerr<< "parsing VCF ... ";
    SnpParser snpFile(params);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // load SV vcf file
    begin = time(NULL);
    std::cerr<< "parsing SV VCF ... ";
    SVParser svFile(params, snpFile);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
 
    //Parse mod vcf file
	begin = time(NULL);
	std::cerr<< "parsing Meth VCF ... ";
    METHParser modFile(params, snpFile, svFile);
	std::cerr<< difftime(time(NULL), begin) << "s\n";
 
    // parsing ref fasta 
    begin = time(NULL);
    std::cerr<< "reading reference ... ";
    std::vector<int> last_pos;
    for(auto chr :snpFile.getChrVec()){
        last_pos.push_back(snpFile.getLastSNP(chr));
    }
    FastaParser fastaParser(params.fastaFile, snpFile.getChrVec(), last_pos, params.numThreads);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // get all detected chromosome
    std::vector<std::string> chrName = snpFile.getChrVec();

    // record all phasing result
    ChrPhasingResult chrPhasingResult;
    // Initialize an empty map in chrPhasingResult to store all phasing results.
    // This is done to prevent issues with multi-threading, by defining an empty map first.
    for (std::vector<std::string>::iterator chrIter = chrName.begin(); chrIter != chrName.end(); chrIter++)    {
        chrPhasingResult[*chrIter] = PhasingResult();
    }

    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};

    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(params.numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
    }

    begin = time(NULL);
    
    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads)
    for(std::vector<std::string>::iterator chrIter = chrName.begin(); chrIter != chrName.end() ; chrIter++ ){
        
        std::time_t chrbegin = time(NULL);
        
        // get last SNP variant position
        int lastSNPpos = snpFile.getLastSNP((*chrIter));
        // therer is no variant on SNP file. 
        if( lastSNPpos == -1 ){
            continue;
        }

	    // fetch chromosome string
        std::string chr_reference = fastaParser.chrString.at(*chrIter);
        // create a bam parser object and prepare to fetch varint from each vcf file
	    BamParser *bamParser = new BamParser((*chrIter), params.bamFile, snpFile, svFile, modFile, chr_reference);
        // use to store variant
        std::vector<ReadVariant> readVariantVec;
        // use to store clip count
        ClipCount clipCount;
        // run fetch variant process
        bamParser->direct_detect_alleles(lastSNPpos, threadPool, params, readVariantVec, clipCount, chr_reference);
        // free memory
        delete bamParser;
        
        // filter variants prone to switch errors in ONT sequencing.
        if(params.isONT){
            snpFile.filterSNP((*chrIter), readVariantVec, chr_reference);
        }

        // bam files are partial file or no read support this chromosome's SNP
        if( readVariantVec.size() == 0 ){
            continue;
        }
        Clip *clip = new Clip((*chrIter), clipCount);
        clip->getCNVInterval(clipCount, (*chrIter));
        

        // create a graph object and prepare to phasing.
        VairiantGraph *vGraph = new VairiantGraph(chr_reference, params, (*chrIter));
        // trans read-snp info to edge info
        vGraph->addEdge(readVariantVec, *clip);
        // run main algorithm
        vGraph->phasingProcess();
        // push result to phasingResult
        vGraph->exportResult((*chrIter), chrPhasingResult[*chrIter]);
        // generate dot file
        if(params.generateDot){
            vGraph->writingDotFile((*chrIter));
        }
        
        // release the memory used by the object.
        vGraph->destroy();
        
        // free memory
        readVariantVec.clear();
        readVariantVec.shrink_to_fit();
        delete vGraph;
        
        std::cerr<< "(" << (*chrIter) << "," << difftime(time(NULL), chrbegin) << "s)";
    }
    hts_tpool_destroy(threadPool.pool);

    std::cerr<< "\nparsing total:  " << difftime(time(NULL), begin) << "s\n";
    
    begin = time(NULL);
    std::cerr<< "merge results ... ";
    // Create a container for merged phasing results.
    PhasingResult mergedPhasingResult;
    // Merge phasing results from all chromosomes.
    mergeAllChrPhasingResult(chrPhasingResult, mergedPhasingResult);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
    
    begin = time(NULL);
    std::cerr<< "writeResult SNP ... ";
    snpFile.writeResult(mergedPhasingResult);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
    
    if(params.svFile!=""){
        begin = time(NULL);
        std::cerr<< "write SV Result ... ";
        svFile.writeResult(mergedPhasingResult);
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }
    
    if(params.modFile!=""){
        begin = time(NULL);
        std::cerr<< "write mod Result ... ";
        modFile.writeResult(mergedPhasingResult);
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }

    std::cerr<< "\ntotal process: " << difftime(time(NULL), processBegin) << "s\n";

    return;
};

PhasingProcess::~PhasingProcess(){
};

