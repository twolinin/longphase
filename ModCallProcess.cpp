#include "ModCallProcess.h"
#include "ModCallParsingBam.h"
#if defined(VALIDATE_PARSECG) || defined(VALIDATE_EXPORT) || \
    defined(VALIDATE_PARSECG_DEPTH) || defined(VALIDATE_JMG)
#include <fstream>
#endif

static MethBamParser* runBamPipeline(
    const std::string &chrName,
    std::string &chrSeq,
    int chrLen,
    ModCallParameters &params,
    htsThreadPool &threadPool,
    MethSnpParser &snpFile,
    std::vector<ModReadVariant> &readVariantVec,
    std::vector<ModReadVariant> &modReadVariantVec)
{
    MethBamParser *methbamparser = new MethBamParser(chrName, params, snpFile, chrSeq);
    // the parser keeps what it needs of the reference
    std::string().swap(chrSeq);
    methbamparser->detectMeth(chrName, chrLen, threadPool, readVariantVec);
#ifdef VALIDATE_PARSECG_DEPTH
    #pragma omp critical (validate_parsecg_depth_dump)
    {
        std::ofstream depthout(params.resultPrefix + "_parsecg_depth.tsv", std::ios::app);
        methbamparser->dumpReadStartEnd(chrName, depthout);
    }
#endif
#ifdef VALIDATE_PARSECG
    #pragma omp critical (validate_parsecg_dump)
    {
        std::ofstream methout(params.resultPrefix + "_parsecg_methmap.tsv", std::ios::app);
        methbamparser->dumpMethMap(chrName, methout);

        std::ofstream readout(params.resultPrefix + "_parsecg_readvar.tsv", std::ios::app);
        for (const auto &rv : readVariantVec) {
            for (const auto &v : rv.variantVec) {
                readout << chrName
                        << "\t" << rv.read_name
                        << "\t" << rv.is_reverse
                        << "\t" << rv.reference_start
                        << "\t" << v.position
                        << "\t" << static_cast<int>(v.allele)
                        << "\t" << v.quality
                        << "\t" << static_cast<int>(v.type)
                        << "\n";
            }
        }
    }
#endif
    methbamparser->calculateDepth();
    methbamparser->judgeMethGenotype(chrName, readVariantVec, modReadVariantVec);
#ifdef VALIDATE_JMG
    #pragma omp critical (validate_jmg_dump)
    {
        std::ofstream genomapout(params.resultPrefix + "_jmg_genomap.tsv", std::ios::app);
        methbamparser->dumpMethGenoMap(chrName, genomapout);
    }
#endif
    return methbamparser;
}

static void runGraphPipeline(
    const std::string &chrName,
    ModCallParameters &params,
    std::vector<ModReadVariant> &modReadVariantVec,
    std::vector<int> &passPosition,
    MethSnpParser &snpFile)
{
    MethylationGraph *modGraph = new MethylationGraph(params);
    modGraph->addEdge(modReadVariantVec, chrName);
    modGraph->connectResults(chrName, passPosition, snpFile.hasValidSnpData());
    modGraph->destroy();
    delete modGraph;
}

static void processChromosome(
    const ReferenceChromosome &chr,
    ModCallParameters &params,
    htsThreadPool &threadPool,
    MethSnpParser &snpFile,
    std::ostringstream &result)
{
    std::time_t chrbegin = time(NULL);

    std::vector<ModReadVariant> readVariantVec;
    std::vector<ModReadVariant> modReadVariantVec;
    std::vector<int> passPosition;

    std::string chrName = chr.name;
    std::string chrSeq = MethFastaParser::fetchSequence(params.fastaFile, chr);
    int chrLen          = chr.length;

    MethBamParser *methbamparser = runBamPipeline(chrName, chrSeq, chrLen, params, threadPool, snpFile, readVariantVec, modReadVariantVec);
    // the graph and the export use only modReadVariantVec and the parser's position table
    std::vector<ModReadVariant>().swap(readVariantVec);
    runGraphPipeline(chrName, params, modReadVariantVec, passPosition, snpFile);
    std::vector<ModReadVariant>().swap(modReadVariantVec);

    methbamparser->exportResult(chrName, chrLen, passPosition, result);
#ifdef VALIDATE_EXPORT
    #pragma omp critical (validate_export_dump)
    {
        std::ofstream out(params.resultPrefix + "_export_dump.tsv", std::ios::app);
        out << result.str();
    }
#endif
    delete methbamparser;

    std::cerr<< "(" << chrName << "," << difftime(time(NULL), chrbegin) << "s)";
}

// jemalloc's mallctl (longphase links jemalloc)
extern "C" int mallctl(const char *name, void *oldp, size_t *oldlenp,
                       void *newp, size_t newlen);

ModCallProcess::ModCallProcess(ModCallParameters params){

    // modcall frees large tables chromosome by chromosome; give freed pages
    // back after 1 s instead of jemalloc's 10 s, so they do not add up to the
    // peak. This lowers modcall's peak by about a quarter at the same speed.
    // Phase is left as it is: there it costs system time for little memory.
    // The default covers the arenas made from now on; the ones in use are set one by one.
    ssize_t decayMs = 1000;
    mallctl("arenas.dirty_decay_ms", nullptr, nullptr, &decayMs, sizeof(decayMs));
    unsigned narenas = 0;
    size_t len = sizeof(narenas);
    if (mallctl("arenas.narenas", &narenas, &len, nullptr, 0) == 0) {
        for (unsigned i = 0; i < narenas; i++) {
            bool initialized = false;
            len = sizeof(initialized);
            std::string arena = "arena." + std::to_string(i);
            if (mallctl((arena + ".initialized").c_str(), &initialized, &len, nullptr, 0) == 0 && initialized)
                mallctl((arena + ".dirty_decay_ms").c_str(), nullptr, nullptr, &decayMs, sizeof(decayMs));
        }
    }

    // load SNP vcf file
    std::time_t begin = time(NULL);
    std::cerr<< "parsing VCF ... ";
    MethSnpParser snpFile(params);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // parsing ref fasta
    begin = time(NULL);
    std::cerr<< "reading reference ... ";
    std::vector<ReferenceChromosome> chrInfo;
    MethFastaParser MethFastaParser(params.fastaFile, chrInfo);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // record all ModCall result; pre-populate to avoid race conditions in parallel region
    std::map<std::string,std::ostringstream> chrModCallResult;
    for (auto chrIter = chrInfo.begin(); chrIter != chrInfo.end(); chrIter++) {
        chrModCallResult[chrIter->name] = std::ostringstream();
    }

    htsThreadPool threadPool = {NULL, 0};
    if (!(threadPool.pool = hts_tpool_init(params.numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
    }

    begin = time(NULL);
    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads)
    for(auto chrIter = chrInfo.begin(); chrIter != chrInfo.end(); ++chrIter) {
        processChromosome(*chrIter, params, threadPool, snpFile, chrModCallResult[chrIter->name]);
    }
    hts_tpool_destroy(threadPool.pool);
    std::cerr<< "\nmodcall total:  " << difftime(time(NULL), begin) << "s\n";

    begin = time(NULL);
    std::cerr<<"write vcf " << " ... ";
    writeResultVCF(params, chrInfo, chrModCallResult);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
}

ModCallProcess::~ModCallProcess(){
};
