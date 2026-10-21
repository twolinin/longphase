#include "PhasingProcess.h"
#include "PhasingGraph.h"
#include "ParsingBam.h"
#include "GNNModel.h"   // compiled-in model weights
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <thread>
#include <set>
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <htslib/hfile.h>
#include <htslib/kseq.h>

// Feature counts come from the generated weight header, so the code and the
// model can never disagree: a mismatch is a compile error.
constexpr int NODE_FEAT_DIM = gnn::kNodeFeat;
constexpr int EDGE_FEAT_DIM = gnn::kEdgeFeat;
constexpr int HP_RADIUS     = 15;
constexpr int GC_RADIUS     = 50;
constexpr int NEIGH_RADIUS  = 500;
constexpr int MAX_NODES     = 256;

PhasingProcess::PhasingProcess(PhasingParameters params)
{
    std::cerr<< "LongPhase Ver " << params.version << "\n";
    std::cerr<< "\n";
    std::cerr<< "--- File Parameter --- \n";
    std::cerr<< "SNP File           : " << params.snpFile      << "\n";
    if (!params.sample.empty())
        std::cerr<< "Sample             : " << params.sample       << "\n";
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
    if (params.phaseIndel) {
        std::cerr<< "Indel Quality      : " << params.indelQuality << "\n";
    }
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
    std::cerr<< "Use GNN            : " << ( params.enableGNN ? "True" : "False" ) << "\n";
    if (params.enableGNN) {
        std::cerr<< "GNN Break Threshold: " << params.gnnBreakThreshold << "\n";
        std::cerr<< "GNN PE Threshold   : " << params.gnnPeThreshold << "\n";
        std::cerr<< "GNN Window         : " << params.gnnWindow << "\n";
        std::cerr<< "GNN Respect Bridge : " << ( params.gnnRespectBridge ? "True" : "False" ) << "\n";
        std::cerr<< "GNN Split Blocks   : " << ( params.gnnSplitBlocks ? "True" : "False" ) << "\n";
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
 
    // get all detected chromosome
    std::vector<std::string> chrName = snpFile.getChrVec();

    // record all phasing result
    ChrPhasingResult chrPhasingResult;
    // Initialize an empty map in chrPhasingResult to store all phasing results.
    // This is done to prevent issues with multi-threading, by defining an empty map first.
    for (std::vector<std::string>::iterator chrIter = chrName.begin(); chrIter != chrName.end(); chrIter++)    {
        chrPhasingResult[*chrIter] = PhasingResult();
    }
    // Graph edges for the GNN, per chromosome; pre-created for the
    // same reason as chrPhasingResult.
    std::map<std::string, GnnEdgeList> chrGnnEdges;
    if (params.enableGNN) {
        for (std::vector<std::string>::iterator chrIter = chrName.begin(); chrIter != chrName.end(); chrIter++) {
            chrGnnEdges[*chrIter];
        }
    }

    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};

    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(params.numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
    }

    begin = time(NULL);
    std::mutex referenceMutex;
    
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

        // The reference is only read near the variants: mark them, then free
        // it. One chromosome's reference is held at a time.
        // (a mutex, not omp critical: waiting threads sleep instead of spinning)
        bool hasReference;
        {
            std::lock_guard<std::mutex> lock(referenceMutex);
            std::string chr_reference = loadChrReference(params.fastaFile, *chrIter, lastSNPpos);
            snpFile.markChrVariants(*chrIter, chr_reference);
            hasReference = !chr_reference.empty();
        }
        // create a bam parser object and prepare to fetch varint from each vcf file
	    BamParser *bamParser = new BamParser((*chrIter), params.bamFile, snpFile, svFile, modFile);
        // use to store variant
        std::vector<ReadVariant> readVariantVec;
        // use to store clip count
        ClipCount clipCount;
        // run fetch variant process
        bamParser->direct_detect_alleles(lastSNPpos, threadPool, params, readVariantVec, clipCount, hasReference);
        // free memory
        delete bamParser;
        
        // filter variants prone to switch errors in ONT sequencing.
        if(params.isONT){
            snpFile.filterSNP((*chrIter), readVariantVec);
        }

        // bam files are partial file or no read support this chromosome's SNP
        if( readVariantVec.size() == 0 ){
            continue;
        }
        Clip *clip = new Clip((*chrIter), clipCount);
        clip->getCNVInterval(clipCount);
        

        // create a graph object and prepare to phasing.
        VairiantGraph *vGraph = new VairiantGraph(params, (*chrIter));
        // trans read-snp info to edge info
        vGraph->addEdge(readVariantVec, *clip);
        // run main algorithm
        vGraph->phasingProcess();
        // push result to phasingResult
        vGraph->exportResult((*chrIter), chrPhasingResult[*chrIter]);
        if(params.enableGNN){
            vGraph->exportGnnEdges(chrGnnEdges.at(*chrIter));
        }
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
        delete clip;

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
    
    // With the GNN these VCFs are intermediate: the GNN reads them and
    // writes the final <prefix>.vcf, then they are removed.
    std::string phasedPrefix = params.enableGNN ? params.resultPrefix + ".pregnn" : params.resultPrefix;

    begin = time(NULL);
    std::cerr<< "writeResult SNP ... ";
    snpFile.writeResult(mergedPhasingResult, phasedPrefix);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
    
    if(params.svFile!=""){
        begin = time(NULL);
        std::cerr<< "write SV Result ... ";
        svFile.writeResult(mergedPhasingResult, phasedPrefix);
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }
    
    if(params.modFile!=""){
        begin = time(NULL);
        std::cerr<< "write mod Result ... ";
        modFile.writeResult(mergedPhasingResult, phasedPrefix);
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }

    if(params.enableGNN){
        // Free the largest phase data before running the GNN.
        PhasingResult().swap(mergedPhasingResult);
        ChrPhasingResult().swap(chrPhasingResult);

        begin = time(NULL);
        std::cerr<< "\nGNN ...\n";
        GNNModule::Params gnnParams;
        gnnParams.vcf_path        = phasedPrefix + ".vcf";
        gnnParams.reference_path  = params.fastaFile;
        gnnParams.sample          = params.sample;
        gnnParams.output_vcf      = params.resultPrefix + ".vcf";
        if(params.svFile!=""){
            gnnParams.sv_vcf         = phasedPrefix + "_SV.vcf";
            gnnParams.output_sv_vcf  = params.resultPrefix + "_SV.vcf";
        }
        if(params.modFile!=""){
            gnnParams.mod_vcf        = phasedPrefix + "_mod.vcf";
            gnnParams.output_mod_vcf = params.resultPrefix + "_mod.vcf";
        }
        gnnParams.break_threshold = params.gnnBreakThreshold;
        gnnParams.pe_threshold    = params.gnnPeThreshold;
        gnnParams.window          = params.gnnWindow;
        gnnParams.threads         = params.numThreads;
        gnnParams.respect_bridge  = params.gnnRespectBridge;
        gnnParams.split_blocks    = params.gnnSplitBlocks;

        GNNModule gnnModule(gnnParams);
        gnnModule.setDotEdges(std::move(chrGnnEdges));
        gnnModule.run();

        std::remove(gnnParams.vcf_path.c_str());
        if(params.svFile!="")  std::remove(gnnParams.sv_vcf.c_str());
        if(params.modFile!="") std::remove(gnnParams.mod_vcf.c_str());
        std::cerr<< "GNN ... " << difftime(time(NULL), begin) << "s\n";
    }

    std::cerr<< "\ntotal process: " << difftime(time(NULL), processBegin) << "s\n";

    return;
};

PhasingProcess::~PhasingProcess(){
};


// ═══════════════════════════════════════════════════════════
// GNN stage of phasing
//
// Optimizations vs the initial version:
//   1. Graph edges are passed from phase in memory (no DOT files)
//   2. Bridge detection: derived per window
//   3. Variant lookup: unordered_map (was linear scan)
//   4. Per-chromosome data precomputed once
// ═══════════════════════════════════════════════════════════

static constexpr float EPS = 1e-6f;

GNNModule::GNNModule(const Params& p)
    : params_(p)
{}
GNNModule::~GNNModule() = default;

// ═══════════════════════════════════════════════════════════
int GNNModule::run() {
    std::cerr << "[GNN] Loading model\n";
    loadModel();
    std::cerr << "[GNN] Parsing VCF\n";
    parseVCF();
    if (!params_.sv_vcf.empty()) {
        std::cerr << "[GNN] Parsing SV VCF\n";
        parseSecondaryVCF(params_.sv_vcf, GVT_SV);
    }
    if (!params_.mod_vcf.empty()) {
        std::cerr << "[GNN] Parsing Methyl VCF\n";
        parseSecondaryVCF(params_.mod_vcf, GVT_METHYL);
    }
    // Re-sort after merging secondary VCFs
    if (!params_.sv_vcf.empty() || !params_.mod_vcf.empty()) {
        for (auto& kv : variants_) {
            auto& vl = kv.second;
            std::sort(vl.begin(), vl.end(),
                      [](const VariantInfo& a, const VariantInfo& b){ return a.pos < b.pos; });
        }
    }
    std::cerr << "[GNN] Loading graph edges\n";
    loadDotEdges();
    if (!params_.reference_path.empty()) {
        std::cerr << "[GNN] Computing genomic features\n";
        computeGenomicFeatures();
    }

    // Pre-compute per-chromosome data
    std::map<std::string, ChromData> chrom_data;
    std::vector<std::string> chroms;
    for (auto& kv : variants_) chroms.push_back(kv.first);

    for (auto& ch : chroms) {
        auto& cd = chrom_data[ch];
        auto& vl = variants_[ch];
        for (auto& v : vl) cd.ps_count[v.ps]++;
        cd.all_pos.reserve(vl.size());
        for (auto& v : vl) cd.all_pos.push_back(v.pos);
        for (auto& v : vl) if (v.indel_len > 0) cd.indel_pos.push_back(v.pos);
        std::sort(cd.indel_pos.begin(), cd.indel_pos.end());

        // Bridge vertices used to be precomputed here over the whole
        // chromosome. They are now derived per window inside
        // processChromosome(), which is what the model was trained on and is
        // also cheaper, since each window holds only a few dozen variants.

        // pos → index in vlist
        for (int i = 0; i < (int)vl.size(); ++i)
            cd.pos_to_idx[vl[i].pos] = i;
    }

    // Build work queue
    struct WorkItem { std::string chrom; int tidx; };
    std::vector<WorkItem> work;
    for (auto& ch : chroms) {
        auto& vl = variants_[ch];
        for (int i = 0; i < (int)vl.size(); ++i)
            if (vl[i].pe >= params_.pe_threshold)
                work.push_back({ch, i});
    }
    std::cerr << "[GNN] " << work.size() << " triggers, " << params_.threads << " threads\n";

    std::atomic<size_t> widx{0}, done{0};
    std::vector<std::thread> threads;
    for (int t = 0; t < params_.threads; ++t) {
        threads.emplace_back([&]() {
            while (true) {
                size_t i = widx.fetch_add(1);
                if (i >= work.size()) return;
                processChromosome(work[i].chrom, work[i].tidx, chrom_data[work[i].chrom]);
                size_t d = done.fetch_add(1) + 1;
                if (d % 2000 == 0)
                    std::cerr << "  " << d << "/" << work.size() << "\r";
            }
        });
    }
    for (auto& t : threads) t.join();
    if (skipped_windows_.load() > 0)
        std::cerr << "\n[GNN] WARNING: " << skipped_windows_.load()
                  << " windows skipped because they exceed " << MAX_NODES
                  << " graph nodes; reduce --window\n";
    shareMergedPredictions();

    int total = 0, unph = 0;
    for (auto& kv_ch : predictions_)
        for (auto& kv_pos : kv_ch.second) {
            auto& pr = kv_pos.second;
            total++;
            if (shouldUnphase(pr)) unph++;
        }
    std::cerr << "\n[GNN] Predictions: " << total << ", unphase: " << unph << "\n";
    if (params_.split_blocks) {
        std::cerr << "[GNN] Checking block connectivity\n";
        computeBlockSplits();
    }
    computeOrphans();
    writeOutputVCF();
    if (!params_.sv_vcf.empty() && !params_.output_sv_vcf.empty())
        writeOutputVCF(params_.sv_vcf, params_.output_sv_vcf);
    if (!params_.mod_vcf.empty() && !params_.output_mod_vcf.empty())
        writeOutputVCF(params_.mod_vcf, params_.output_mod_vcf);
    return 0;
}

// ═══════════════════════════════════════════════════════════
void GNNModule::loadModel() {
    // Decodes the weights compiled in via GNNWeights.h. Nothing is read from
    // disk, so --model is not consulted.
    model_.reset(new gnn::Model());
    std::cerr << "[GNN]   " << gnn::kParamCount << " parameters, "
              << NODE_FEAT_DIM << " node / " << EDGE_FEAT_DIM
              << " edge features\n";
}

// ═══════════════════════════════════════════════════════════
// The sample of the SNP VCF to read and write: --sample, by default the
// first one (phase has checked that the sample exists).
int GNNModule::sampleIndex(bcf_hdr_t* hdr) const {
    if (!params_.sample.empty()) {
        int i = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, params_.sample.c_str());
        if (i >= 0) return i;
    }
    return 0;
}

// Set the PS of sample si; missing removes it. With one sample the PS field
// is removed from the record, as before; with several, only this sample's
// value changes.
void GNNModule::setSamplePS(bcf_hdr_t* hdr, bcf1_t* rec, int si, int32_t ps) {
    const int nsmpl = bcf_hdr_nsamples(hdr);
    if (nsmpl == 1) {
        if (ps == bcf_int32_missing) bcf_update_format_int32(hdr, rec, "PS", nullptr, 0);
        else bcf_update_format_int32(hdr, rec, "PS", &ps, 1);
        return;
    }
    int32_t* v = nullptr; int nv = 0;
    if (bcf_get_format_int32(hdr, rec, "PS", &v, &nv) != nsmpl) {
        free(v);
        v = (int32_t*)malloc(nsmpl * sizeof(int32_t));
        for (int i = 0; i < nsmpl; ++i) v[i] = bcf_int32_missing;
    }
    v[si] = ps;
    bcf_update_format_int32(hdr, rec, "PS", v, nsmpl);
    free(v);
}

namespace {
// readVcfRecords: records per batch, and threads; reading and writing the
// records on one thread limit the gain from more threads
constexpr size_t kVcfBatch = 1024;
constexpr int kVcfThreads = 8;
}

void GNNModule::parseVCF() {
    htsFile* fp = hts_open(params_.vcf_path.c_str(), "r");
    if (!fp) {
        std::cerr << "[GNN] ERROR: cannot open " << params_.vcf_path << "\n";
        std::exit(EXIT_FAILURE);
    }
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        std::cerr << "[GNN] ERROR: cannot read VCF header of " << params_.vcf_path << "\n";
        std::exit(EXIT_FAILURE);
    }
    int n = 0;
    const int si = sampleIndex(hdr);
    // the variant of each record of a batch, if it is a phased heterozygote
    std::vector<VariantInfo> got(kVcfBatch);
    std::vector<char> keep(kVcfBatch);
    readVcfRecords(fp, hdr, std::min(params_.threads, kVcfThreads), kVcfBatch,
        [&](size_t i, bcf_hdr_t* hdr, bcf1_t* rec) {
        keep[i] = 0;
        bcf_unpack(rec, BCF_UN_ALL);
        int32_t* gt = nullptr; int ng = 0;
        if (bcf_get_genotypes(hdr, rec, &gt, &ng) < 0) return;
        const int ploidy = ng / bcf_hdr_nsamples(hdr);
        if (ploidy < 2) { free(gt); return; }
        int32_t* sg = gt + si * ploidy;
        if (!bcf_gt_is_phased(sg[1])) { free(gt); return; }
        int a0 = bcf_gt_allele(sg[0]), a1 = bcf_gt_allele(sg[1]);
        free(gt);
        if (a0 == a1) return;

        VariantInfo& vi = got[i];
        vi = VariantInfo();
        vi.chrom = bcf_hdr_id2name(hdr, rec->rid);
        vi.pos = rec->pos; vi.gt_ref = a0; vi.gt_alt = a1;
        vi.is_phased = true; vi.pe = vi.h1 = vi.h2 = 0; vi.ps = -1;
        int rl = (int)strlen(rec->d.allele[0]);
        int al = rec->n_allele > 1 ? (int)strlen(rec->d.allele[1]) : rl;
        vi.indel_len = abs(al - rl);

        float* fv = nullptr; int nv = 0;
        if (bcf_get_info_float(hdr, rec, "PE", &fv, &nv) > 0) { vi.pe = fv[0]; free(fv); fv=nullptr; nv=0; }
        if (bcf_get_info_float(hdr, rec, "H1", &fv, &nv) > 0) { vi.h1 = fv[0]; free(fv); fv=nullptr; nv=0; }
        if (bcf_get_info_float(hdr, rec, "H2", &fv, &nv) > 0) { vi.h2 = fv[0]; free(fv); fv=nullptr; nv=0; }
        int32_t* ps = nullptr; int np = 0;
        if (bcf_get_format_int32(hdr, rec, "PS", &ps, &np) > si) vi.ps = ps[si];
        free(ps);
        vi.hp_len = vi.gc_content = vi.str_context = vi.seq_entropy = 0;
        vi.vtype = (vi.indel_len == 0) ? GVT_SNV : GVT_INDEL;
        keep[i] = 1;
    },
        [&](size_t i, bcf_hdr_t*, bcf1_t*) {
        if (!keep[i]) return;
        variants_[got[i].chrom].push_back(got[i]);
        n++;
    });
    bcf_hdr_destroy(hdr); hts_close(fp);
    for (auto& kv : variants_) {
        auto& vl = kv.second;
        std::sort(vl.begin(), vl.end(),
                  [](const VariantInfo& a, const VariantInfo& b){ return a.pos < b.pos; });
    }
    std::cerr << "  " << n << " phased variants\n";
}

// ═══════════════════════════════════════════════════════════
// Graph edges from phase
// ═══════════════════════════════════════════════════════════
void GNNModule::setDotEdges(std::map<std::string, GnnEdgeList>&& edges) {
    pending_edges_ = std::move(edges);
}

void GNNModule::loadDotEdges() {
    // Take over the edges of each chromosome with variants
    size_t total = 0;
    for (auto& kv : variants_) {
        const std::string& ch = kv.first;
        auto& em = dot_edges_[ch];
        auto pit = pending_edges_.find(ch);
        if (pit == pending_edges_.end()) continue;
        em.swap(pit->second);
        em.finish();
        total += em.edgeCount();
    }
    pending_edges_.clear();
    std::cerr << "  " << total << " graph edges\n";
}

void GnnEdgeList::finish() {
    if (sorted) return;
    // order the groups by source position; a stable sort keeps the groups
    // of one source in the order they were added
    std::vector<size_t> order(src.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(),
                     [&](size_t x, size_t y) { return src[x] < src[y]; });
    GnnEdgeList out;
    out.links.reserve(links.size());
    for (size_t i : order) {
        size_t b = start[i], e = i + 1 < start.size() ? start[i + 1] : links.size();
        out.src.push_back(src[i]);
        out.weight.push_back(weight[i]);
        out.start.push_back(out.links.size());
        out.links.insert(out.links.end(), links.begin() + b, links.begin() + e);
    }
    swap(out);
    sorted = true;
}

// ═══════════════════════════════════════════════════════════
void GNNModule::computeGenomicFeatures() {
    std::vector<std::string> chroms;
    for (auto& kv : variants_) chroms.push_back(kv.first);

    std::atomic<size_t> cidx{0};
    int n_threads = std::min(params_.threads, (int)chroms.size());
    std::vector<std::thread> threads;

    for (int t = 0; t < n_threads; ++t) {
        threads.emplace_back([&]() {
            // Each thread gets its own reader
            FastaReader reader(params_.reference_path);
            if (!reader.loaded()) return;
            while (true) {
                size_t i = cidx.fetch_add(1);
                if (i >= chroms.size()) return;
                auto& ch = chroms[i];
                for (auto& v : variants_[ch]) {
                    // The model was trained with these four features set
                    // to 0 for SVs, whose breakpoint context is not
                    // meaningful, so leave them at 0 here too.
                    if (v.vtype == GVT_SV) continue;
                    int len = 0;
                    int s1 = std::max(0, v.pos - HP_RADIUS);
                    std::string ss = reader.fetch(ch, s1, v.pos+HP_RADIUS, len, true);
                    if (len > 0) {
                        v.seq_entropy = calcSeqEntropy(ss);
                        v.hp_len = (float)calcHPLength(ss, v.pos-s1);
                        v.str_context = calcSTRContext(ss, v.pos-s1);
                    }
                    int s2 = std::max(0, v.pos - GC_RADIUS);
                    ss = reader.fetch(ch, s2, v.pos+GC_RADIUS, len, true);
                    if (len > 0) {
                        v.gc_content = calcGCContent(ss);
                    }
                }
            }
        });
    }
    for (auto& t : threads) t.join();
}

// ═══════════════════════════════════════════════════════════
// Process one trigger
// ═══════════════════════════════════════════════════════════
void GNNModule::processChromosome(const std::string& chrom, int tidx,
                                  const ChromData& cd) {
    auto& vlist = variants_[chrom];
    auto& edge_map = dot_edges_[chrom];
    auto& ps_count = cd.ps_count;
    auto& all_pos = cd.all_pos;
    auto& indel_pos = cd.indel_pos;

    int center_pos = vlist[tidx].pos;
    int center_ps  = vlist[tidx].ps;
    int wi_start = std::max(0, tidx - params_.window);
    int wi_end   = std::min((int)vlist.size() - 1, tidx + params_.window);

    std::vector<const VariantInfo*> wvars;
    // the window's positions, sorted and distinct; a position's index in wpos
    // keys the per-position values below
    std::vector<int> wpos;
    for (int i = wi_start; i <= wi_end; ++i) {
        wvars.push_back(&vlist[i]);
        wpos.push_back(vlist[i].pos);
    }
    std::sort(wpos.begin(), wpos.end());
    wpos.erase(std::unique(wpos.begin(), wpos.end()), wpos.end());
    const int nk = (int)wpos.size();
    // index of a position in wpos, -1 outside the window
    auto wk = [&](int pos) -> int {
        auto it = std::lower_bound(wpos.begin(), wpos.end(), pos);
        return (it != wpos.end() && *it == pos) ? (int)(it - wpos.begin()) : -1;
    };
    int n_var = (int)wvars.size();
    int N = n_var * 2;
    if (N > MAX_NODES) { skipped_windows_.fetch_add(1); return; }
    // the position index of each window variant; the last variant at a
    // position stands for it
    std::vector<int> kv(n_var);
    std::vector<const VariantInfo*> pos_to_var(nk);
    for (int i = 0; i < n_var; ++i) {
        kv[i] = wk(wvars[i]->pos);
        pos_to_var[kv[i]] = wvars[i];
    }
    // the edges of each position
    std::vector<char> has_edges(nk);
    std::vector<GnnEdgeList::Range> edges_of(nk);
    for (int k = 0; k < nk; ++k) {
        has_edges[k] = edge_map.count(wpos[k]) > 0;
        if (has_edges[k]) edges_of[k] = edge_map.at(wpos[k]);
    }

    float max_offset = 1.0f;
    for (auto* v : wvars) {
        float d = std::abs((float)(v->pos - center_pos));
        if (d > max_offset) max_offset = d;
    }

    // h1/h2 from DOT first edge per allele
    std::vector<char> has_h1(nk, 0), has_h2(nk, 0);
    std::vector<float> dot_h1(nk), dot_h2(nk);
    for (int i = 0; i < n_var; ++i) {
        const int k = kv[i];
        if (!has_edges[k]) continue;
        for (const GnnEdge& e : edges_of[k]) {
            if (e.src_allele == 1 && !has_h1[k]) { dot_h1[k] = e.weight; has_h1[k] = 1; }
            else if (e.src_allele == 2 && !has_h2[k]) { dot_h2[k] = e.weight; has_h2[k] = 1; }
        }
    }

    // out_degree, out_weight_std
    std::vector<float> out_deg(nk), out_wstd(nk);
    std::vector<int> dsts; std::vector<float> ws;
    for (int i = 0; i < n_var; ++i) {
        const int k = kv[i];
        if (!has_edges[k]) { out_deg[k]=0; out_wstd[k]=0; continue; }
        dsts.clear(); ws.clear();
        for (const GnnEdge& e : edges_of[k]) { dsts.push_back(e.dst_pos); ws.push_back(e.weight); }
        // the number of distinct targets
        std::sort(dsts.begin(), dsts.end());
        const size_t n_dst = std::unique(dsts.begin(), dsts.end()) - dsts.begin();
        out_deg[k] = std::min((float)n_dst/10.0f, 1.0f);
        if (ws.size() > 1) {
            float mu=0; for(float w:ws) mu+=w; mu/=ws.size();
            float var=0; for(float w:ws) var+=(w-mu)*(w-mu); var/=ws.size();
            out_wstd[k] = std::min(std::sqrt(var)/100.0f, 1.0f);
        } else out_wstd[k] = 0;
    }

    // window_mean_votes (VCF h1+h2, full chromosome ±10)
    std::vector<float> wmean(nk);
    for (int i = wi_start; i <= wi_end; ++i) {
        int lo = std::max(0, i-10), hi = std::min((int)vlist.size()-1, i+10);
        float sum=0; int cnt=0;
        for (int j=lo; j<=hi; ++j) { sum += vlist[j].h1 + vlist[j].h2; cnt++; }
        wmean[kv[i - wi_start]] = cnt > 0 ? sum/cnt : 1.0f;
    }

    // Neighbor density lambda
    auto neigh_dens = [&](int pos) -> float {
        auto lo = std::lower_bound(all_pos.begin(), all_pos.end(), pos - NEIGH_RADIUS);
        auto hi = std::upper_bound(all_pos.begin(), all_pos.end(), pos + NEIGH_RADIUS);
        return ((float)(hi-lo) - 1.0f) / 20.0f;
    };

    // Nearest indel lambda
    auto dist_indel = [&](int pos) -> float {
        if (indel_pos.empty()) return std::log10(10001.0f);
        auto it = std::lower_bound(indel_pos.begin(), indel_pos.end(), pos);
        float best = 10000.0f;
        if (it != indel_pos.end()) best = std::min(best, (float)std::abs(*it - pos));
        if (it != indel_pos.begin()) { --it; best = std::min(best, (float)std::abs(*it - pos)); }
        return std::min(std::log10(1.0f+best), std::log10(10001.0f));
    };

    // ── Bridge vertices ──────────────────────────────────
    // Mirrors is_bridge_vertex() in prepare_gnn_data.py: the graph is built
    // from edges whose endpoints both lie inside this window, collapsed to
    // position level and treated as undirected. A variant is a bridge when
    // removing it raises the number of connected components of the window.
    //
    // The window holds a few dozen variants, so recomputing per window is
    // cheaper than the chromosome-wide pass this replaces, and it is what
    // the model was trained on.
    // pos_adj[a * nk + b]: positions a and b are joined; dst_k: the target's
    // position index of each edge, in edge order, -1 outside the window
    std::vector<char> pos_adj((size_t)nk * nk, 0), has_nb(nk, 0);
    std::vector<std::vector<int>> dst_k(nk);
    for (int a = 0; a < nk; ++a) {
        if (!has_edges[a]) continue;
        for (const GnnEdge& e : edges_of[a]) {
            const int b = wk(e.dst_pos);
            dst_k[a].push_back(b);
            if (b < 0 || b == a) continue;
            pos_adj[(size_t)a * nk + b] = pos_adj[(size_t)b * nk + a] = 1;
            has_nb[a] = has_nb[b] = 1;
        }
    }

    // Connected components of the window, optionally with one position removed.
    std::vector<char> vis(nk);
    std::vector<int> stack;
    auto n_components = [&](int exclude) -> int {
        std::fill(vis.begin(), vis.end(), 0);
        int count = 0;
        for (int start = 0; start < nk; ++start) {
            if (start == exclude || vis[start]) continue;
            ++count;
            stack.assign(1, start);
            while (!stack.empty()) {
                int node = stack.back(); stack.pop_back();
                if (vis[node]) continue;
                vis[node] = 1;
                for (int nb = 0; nb < nk; ++nb)
                    if (pos_adj[(size_t)node * nk + nb] && nb != exclude && !vis[nb])
                        stack.push_back(nb);
            }
        }
        return count;
    };

    const int comp_before = n_components(-1);
    std::vector<char> bridge(nk, 0);
    for (int k = 0; k < nk; ++k) {
        if (nk <= 1) break;
        // Needs at least one neighbour that survives the removal.
        if (!has_nb[k]) continue;
        if (n_components(k) > comp_before) bridge[k] = 1;
    }

    // out_wsum for weight_ratio, per position and source allele
    std::vector<std::array<float, 2>> out_wsum(nk, std::array<float, 2>{{0.0f, 0.0f}});
    for (int a = 0; a < nk; ++a) {
        if (!has_edges[a]) continue;
        size_t t = 0;
        for (const GnnEdge& e : edges_of[a]) {
            if (dst_k[a][t++] < 0) continue;
            out_wsum[a][e.src_allele - 1] += e.weight;
        }
    }

    // ── Node features [N × NODE_FEAT_DIM] ──
    std::vector<float> nf(N * NODE_FEAT_DIM, 0.0f);
    // node of each position and allele: the last variant at the position
    std::vector<std::array<int, 2>> pa_to_nid(nk, std::array<int, 2>{{-1, -1}});

    for (int i = 0; i < n_var; ++i) {
        const auto* v = wvars[i];
        const int k = kv[i];
        float h1v = has_h1[k] ? dot_h1[k] : v->h1;
        float h2v = has_h2[k] ? dot_h2[k] : v->h2;
        float total = h1v + h2v;
        float vote_r = total > 0 ? std::max(h1v,h2v)/total : 0.0f;
        float log_tot = std::log2(1.0f+total);
        float lcd = std::log10(1.0f + std::abs((float)(v->pos-center_pos)));
        float ps_blk = std::log1p((float)(ps_count.count(v->ps)?ps_count.at(v->ps):1))
                       / std::log1p(100.0f);
        float wm = wmean[k];
        float rvd = std::log2(total / (wm+EPS) + EPS);
        float bridge_f = bridge[k] ? 1.0f : 0.0f;

        for (int a = 0; a < 2; ++a) {
            int allele = a+1, nid = i*2+a;
            pa_to_nid[k][a] = nid;
            int hp_val = (allele==1) ? v->gt_ref : v->gt_alt;
            float* f = &nf[nid * NODE_FEAT_DIM];
            f[0]=v->pe; f[1]=(float)hp_val; f[2]=(float)v->gt_ref; f[3]=(float)v->gt_alt;
            f[4]=(v->pos==center_pos)?1.0f:0.0f;
            f[5]=(v->ps==center_ps && v->ps!=-1)?1.0f:0.0f;
            f[6]=(float)a; f[7]=(float)(v->pos-center_pos)/max_offset;
            f[8]=h1v; f[9]=h2v; f[10]=total; f[11]=vote_r;
            f[12]=log_tot; f[13]=std::log1p(std::max(0.0f,h1v)); f[14]=std::log1p(std::max(0.0f,h2v));
            f[15]=v->hp_len; f[16]=v->gc_content; f[17]=v->str_context; f[18]=v->seq_entropy;
            f[19]=neigh_dens(v->pos); f[20]=out_deg[k]; f[21]=ps_blk;
            f[22]=lcd; f[23]=dist_indel(v->pos); f[24]=rvd;
            f[25]=out_wstd[k]; f[26]=bridge_f;
            // v21 cophasing: variant type indicators
            f[27]=(v->vtype==GVT_SNV)?1.0f:0.0f;
            f[28]=(v->vtype==GVT_INDEL)?1.0f:0.0f;
            f[29]=(v->vtype==GVT_SV)?1.0f:0.0f;
            f[30]=(v->vtype==GVT_METHYL)?1.0f:0.0f;
        }
    }

    // ── Edges ──
    std::vector<float> adj(N*N, 0.0f);
    std::vector<float> ef(N*N*EDGE_FEAT_DIM, 0.0f);

    // A cell (row, col) is written only by edges of one source position, in
    // their order, so the order of the positions does not matter.
    for (int a = 0; a < nk; ++a) {
        if (!has_edges[a]) continue;
        const int pos = wpos[a];
        size_t t = 0;
        for (const GnnEdge& e : edges_of[a]) {
            const int b = dst_k[a][t++];
            if (b < 0) continue;
            int sn = pa_to_nid[a][e.src_allele - 1], dn = pa_to_nid[b][e.dst_allele - 1];
            if (sn < 0 || dn < 0) continue;
            float same_hp = (e.src_allele == e.dst_allele) ? 1.0f : 0.0f;
            auto* sv = pos_to_var[a]; auto* dv = pos_to_var[b];
            float inter_ps = (sv && dv && sv->ps != dv->ps) ? 1.0f : 0.0f;
            float wsum = out_wsum[a][e.src_allele - 1];
            int row = dn, col = sn;
            adj[row*N+col] = 1.0f;
            float* ep = &ef[(row*N+col)*EDGE_FEAT_DIM];
            ep[0] = 1.0f/(1.0f+std::exp(-e.weight));
            ep[1] = std::log10(1.0f+std::abs((float)(e.dst_pos-pos)));
            ep[2] = std::log1p(e.weight);
            ep[3] = same_hp; ep[4] = inter_ps;
            ep[5] = e.weight/(wsum+EPS);
        }
    }

    // Undirected
    for (int i=0;i<N;++i) for (int j=i+1;j<N;++j) {
        bool ij=adj[i*N+j]>0.5f, ji=adj[j*N+i]>0.5f;
        if (ij&&!ji) { adj[j*N+i]=1; for(int f=0;f<EDGE_FEAT_DIM;++f) ef[(j*N+i)*EDGE_FEAT_DIM+f]=ef[(i*N+j)*EDGE_FEAT_DIM+f]; }
        else if (ji&&!ij) { adj[i*N+j]=1; for(int f=0;f<EDGE_FEAT_DIM;++f) ef[(i*N+j)*EDGE_FEAT_DIM+f]=ef[(j*N+i)*EDGE_FEAT_DIM+f]; }
    }

    // Self-loops (per-node mean)
    for (int i=0;i<N;++i) {
        adj[i*N+i]=1.0f;
        float se[EDGE_FEAT_DIM]={}; int ni=0;
        for (int j=0;j<N;++j) if (j!=i && adj[i*N+j]>0.5f) {
            for(int f=0;f<EDGE_FEAT_DIM;++f) se[f]+=ef[(i*N+j)*EDGE_FEAT_DIM+f];
            ni++;
        }
        if (ni>0) for(int f=0;f<EDGE_FEAT_DIM;++f) ef[(i*N+i)*EDGE_FEAT_DIM+f]=se[f]/ni;
    }

    // ── Inference ──
    // adj(row, col) != 0 means an edge col -> row, which is the orientation
    // gnn::Model expects.
    std::vector<float> probs = model_->forward(nf, adj, ef, N);
    const float* op = probs.data();

    // ── Collect predictions for center-zone nodes ──
    // Match Python _center_mask: thr = min(10 / max(1, round(max_r * 20)), 1.0)
    float max_r = 0;
    for (int i = 0; i < n_var; ++i) {
        float r = std::abs((float)(wvars[i]->pos - center_pos) / max_offset);
        if (r > max_r) max_r = r;
    }
    if (max_r < 1e-9f) max_r = 1.0f;
    float center_thr = std::min(10.0f / std::max(1.0f, std::round(max_r * 20.0f)), 1.0f);

    for (int i=0;i<n_var;++i) {
        float rel = std::abs((float)(wvars[i]->pos - center_pos) / max_offset);
        if (rel > center_thr) continue;
        int n1=i*2, n2=i*2+1;
        float p_err = (op[n1*2+1]+op[n2*2+1])/2.0f;
        bool br = bridge[kv[i]] != 0;
        std::lock_guard<std::mutex> lk(pred_mutex_);
        auto& pr = predictions_[chrom][wvars[i]->pos];
        pr.prob_error = (pr.prob_error*pr.count + p_err)/(pr.count+1);
        pr.is_bridge = pr.is_bridge || br;
        pr.count++;
    }
}

// ═══════════════════════════════════════════════════════════
// Block connectivity check: after unphasing, if a PS block is split
// into disconnected components (bridge removed), assign new PS IDs.
// Parallelized per chromosome.
void GNNModule::computeBlockSplits() {
    std::vector<std::string> chroms;
    for (auto& kv : variants_) chroms.push_back(kv.first);

    std::atomic<size_t> cidx{0};
    std::mutex reassign_mutex;
    std::atomic<int> total_splits{0}, total_new{0}, total_reassigned{0};

    auto worker = [&]() {
        while (true) {
            size_t ci = cidx.fetch_add(1);
            if (ci >= chroms.size()) return;
            const std::string& ch = chroms[ci];
            auto& vl = variants_[ch];
            auto& em = dot_edges_[ch];

            // Which positions are being unphased?
            std::unordered_set<int> unphase_pos;
            if (predictions_.count(ch)) {
                for (auto& kv_pos : predictions_.at(ch)) {
                    const auto& pos = kv_pos.first;
                    auto& pr = kv_pos.second;
                    if (shouldUnphase(pr))
                        unphase_pos.insert(pos);
                }
            }
            if (unphase_pos.empty()) continue;

            // Group remaining phased variants by PS
            std::unordered_map<int, std::vector<int>> ps_blocks; // ps -> positions
            std::unordered_map<int, int> pos_to_ps;
            for (auto& v : vl) {
                if (!v.is_phased) continue;
                if (v.ps < 0) continue;
                if (unphase_pos.count(v.pos)) continue;
                ps_blocks[v.ps].push_back(v.pos);
                pos_to_ps[v.pos] = v.ps;
            }

            std::unordered_map<int, int> local_reassign;
            int local_splits = 0, local_new = 0;

            // components: positions of each connected part of a PS block, in
            // the order the parts are first reached from the block's positions
            auto split = [&](std::vector<std::vector<int>>& components) {
                if (components.size() <= 1) return;
                // Sort: largest component keeps original PS
                std::sort(components.begin(), components.end(),
                          [](const std::vector<int>& a, const std::vector<int>& b){
                              return a.size() > b.size(); });
                local_splits++;
                for (size_t i = 1; i < components.size(); ++i) {
                    // New PS = min position in component + 1 (1-based convention)
                    int new_ps = *std::min_element(components[i].begin(),
                                                   components[i].end()) + 1;
                    local_new++;
                    for (int pos : components[i])
                        local_reassign[pos] = new_ps;
                }
            };

            // Graph nodes are positions, joined by same-PS, both-phased edges
            // (and 5mC links); a position's PS is that of its last record.
            std::vector<int> nodes;
            nodes.reserve(pos_to_ps.size());
            for (auto& kv : pos_to_ps) nodes.push_back(kv.first);
            std::sort(nodes.begin(), nodes.end());
            auto node_of = [&](int pos) {
                auto it = std::lower_bound(nodes.begin(), nodes.end(), pos);
                return (it != nodes.end() && *it == pos) ? (int)(it - nodes.begin()) : -1;
            };
            std::vector<int> node_ps(nodes.size());
            for (size_t i = 0; i < nodes.size(); ++i) node_ps[i] = pos_to_ps[nodes[i]];
            // A position with records in several blocks is "stale" in the blocks
            // other than its PS: there it is connected only through edges to
            // other positions that are stale in the same block.
            std::vector<char> stale(nodes.size(), 0);
            for (auto& kv_ps : ps_blocks)
                for (int pos : kv_ps.second) {
                    const int i = node_of(pos);
                    if (node_ps[i] != kv_ps.first) stale[i] = 1;
                }

            std::vector<int> parent(nodes.size());
            for (size_t i = 0; i < parent.size(); ++i) parent[i] = (int)i;
            auto find = [](std::vector<int>& par, int i) {
                while (par[i] != i) i = par[i] = par[par[i]];
                return i;
            };
            auto unite = [&](std::vector<int>& par, int i, int j) {
                i = find(par, i); j = find(par, j);
                if (i != j) par[std::max(i, j)] = std::min(i, j);
            };
            std::vector<std::pair<int, int>> stale_edges;   // both ways
            auto join = [&](int s, int d) {
                if (s < 0 || d < 0 || node_ps[s] != node_ps[d]) return;
                unite(parent, s, d);
                if (stale[s] && stale[d]) {
                    stale_edges.emplace_back(s, d);
                    stale_edges.emplace_back(d, s);
                }
            };
            for (size_t gi = 0; gi < em.groups(); ++gi) {
                const int s = node_of(em.src[gi]);
                if (s < 0) continue;
                for (const GnnEdge& e : em.group(gi)) join(s, node_of(e.dst_pos));
            }
            for (auto& v : vl)
                if (v.link_pos >= 0) join(node_of(v.pos), node_of(v.link_pos));
            std::sort(stale_edges.begin(), stale_edges.end());

            // Per-PS connected components, numbered in the order the block's
            // positions first reach them
            std::vector<int> comp_of(nodes.size(), -1);  // root -> component (roots differ between blocks)
            std::vector<int> seen(nodes.size(), -1);     // last block that listed the position
            std::vector<int> loc(nodes.size(), -1);      // index among the block's stale positions
            int bi = 0;
            for (auto& kv_ps : ps_blocks) {
                auto& positions = kv_ps.second;
                ++bi;
                if (positions.size() <= 1) continue;
                std::vector<int> order, stale_here;
                for (int pos : positions) {
                    const int i = node_of(pos);
                    if (seen[i] == bi) continue;
                    seen[i] = bi;
                    order.push_back(i);
                    if (node_ps[i] != kv_ps.first) { loc[i] = (int)stale_here.size(); stale_here.push_back(i); }
                }
                std::vector<int> lpar(stale_here.size()), lcomp(stale_here.size(), -1);
                for (size_t k = 0; k < lpar.size(); ++k) lpar[k] = (int)k;
                for (int i : stale_here)
                    for (auto it = std::lower_bound(stale_edges.begin(), stale_edges.end(), std::make_pair(i, INT_MIN));
                         it != stale_edges.end() && it->first == i; ++it)
                        if (loc[it->second] >= 0) unite(lpar, loc[i], loc[it->second]);
                std::vector<std::vector<int>> components;
                for (int i : order) {
                    int& c = loc[i] < 0 ? comp_of[find(parent, i)] : lcomp[find(lpar, loc[i])];
                    if (c < 0) { c = (int)components.size(); components.emplace_back(); }
                    components[c].push_back(nodes[i]);
                }
                for (int i : stale_here) loc[i] = -1;
                split(components);
            }

            if (!local_reassign.empty()) {
                std::lock_guard<std::mutex> lk(reassign_mutex);
                ps_reassign_[ch] = std::move(local_reassign);
                total_splits += local_splits;
                total_new += local_new;
                total_reassigned += (int)ps_reassign_[ch].size();
            }
        }
    };

    std::vector<std::thread> threads;
    for (int t = 0; t < params_.threads; ++t) threads.emplace_back(worker);
    for (auto& t : threads) t.join();

    std::cerr << "[GNN] PS blocks split: " << total_splits.load()
              << ", new blocks: " << total_new.load()
              << ", variants reassigned: " << total_reassigned.load() << "\n";
}

// phase merges a run of consecutive 5mC positions into one graph node, so
// the sites of a run were phased as one and must be unphased as one. The
// node's own position carries the graph edges, so its prediction decides
// for every site of the run.
void GNNModule::shareMergedPredictions() {
    int n = 0;
    for (auto& kv : variants_) {
        const std::string& ch = kv.first;
        auto pit = predictions_.find(ch);
        if (pit == predictions_.end()) continue;
        auto& preds = pit->second;
        std::unordered_set<int> phased_pos;
        for (auto& v : kv.second) phased_pos.insert(v.pos);
        for (auto& v : kv.second) {
            if (v.link_pos < 0 || !phased_pos.count(v.link_pos)) continue;
            auto it = preds.find(v.link_pos);
            if (it != preds.end()) {
                Prediction shared = it->second;   // copy before a possible rehash
                preds[v.pos] = shared;
            }
            else preds.erase(v.pos);
            n++;
        }
    }
    std::cerr << "[GNN] 5mC sites sharing their merged node's prediction: " << n << "\n";
}

bool GNNModule::shouldUnphase(const Prediction& pr) const {
    return pr.prob_error >= params_.break_threshold &&
           (!params_.respect_bridge || !pr.is_bridge);
}

// A phase set holding one variant carries no phase information. When GNN
// correction (unphasing, and block splitting if enabled) leaves a variant
// alone in its phase set, unphase it too. Phase sets that already held a
// single variant before correction are left as phase wrote them. SNVs,
// indels, SVs and 5mC sites share phase sets, so they are counted together.
void GNNModule::computeOrphans() {
    int total = 0;
    for (auto& kv : variants_) {
        const std::string& ch = kv.first;
        auto& vl = kv.second;
        auto pit = predictions_.find(ch);
        auto rit = ps_reassign_.find(ch);

        std::unordered_map<int, int> before, after;   // PS -> variants
        struct Kept { int pos, ps_before, ps_after; };
        std::vector<Kept> kept;
        for (auto& v : vl) {
            if (!v.is_phased || v.ps < 0) continue;
            before[v.ps]++;
            if (pit != predictions_.end()) {
                auto it = pit->second.find(v.pos);
                if (it != pit->second.end() && shouldUnphase(it->second)) continue;
            }
            int ps = v.ps;
            if (rit != ps_reassign_.end()) {
                auto it = rit->second.find(v.pos);
                if (it != rit->second.end()) ps = it->second;
            }
            after[ps]++;
            kept.push_back({v.pos, v.ps, ps});
        }
        for (auto& k : kept) {
            if (after[k.ps_after] != 1) continue;
            // Unchanged single-variant phase set from phase: keep it.
            if (k.ps_after == k.ps_before && before[k.ps_before] == 1) continue;
            orphans_[ch].insert(k.pos);
            total++;
        }
    }
    std::cerr << "[GNN] Variants left alone in their phase set: " << total << "\n";
}



// ═══════════════════════════════════════════════════════════
void GNNModule::writeOutputVCF() {
    writeOutputVCF(params_.vcf_path, params_.output_vcf);
}

void GNNModule::writeOutputVCF(const std::string& in_path, const std::string& out_path) {
    htsFile* ifp = hts_open(in_path.c_str(), "r");
    if (!ifp) {
        std::cerr << "[GNN] ERROR: cannot open " << in_path << "\n";
        std::exit(EXIT_FAILURE);
    }
    bcf_hdr_t* hdr = bcf_hdr_read(ifp);
    if (!hdr) {
        std::cerr << "[GNN] ERROR: cannot read VCF header of " << in_path << "\n";
        std::exit(EXIT_FAILURE);
    }
    htsFile* ofp = hts_open(out_path.c_str(), "w");
    if (!ofp) {
        std::cerr << "[GNN] ERROR: cannot create " << out_path << "\n";
        std::exit(EXIT_FAILURE);
    }
    if (bcf_hdr_write(ofp, hdr) < 0) {
        std::cerr << "[GNN] ERROR: failed to write VCF header to " << out_path << "\n";
        std::exit(EXIT_FAILURE);
    }
    int nu = 0, nsplit = 0, norphan = 0;
    // --sample applies to the SNP VCF; the SV and 5mC VCFs written by phase
    // carry the phase in their first sample column
    const int si = (in_path == params_.vcf_path) ? sampleIndex(hdr) : 0;
    const int nsmpl = bcf_hdr_nsamples(hdr);
    // Each record of a batch is changed and formatted as VCF text (as
    // bcf_write does) on the threads, then written in order.
    struct Out {
        kstring_t text = KS_INITIALIZE;
        bool formatted;
        char unph, orphan, split;
    };
    std::vector<Out> out(kVcfBatch);
    readVcfRecords(ifp, hdr, std::min(params_.threads, kVcfThreads), kVcfBatch,
        [&](size_t i, bcf_hdr_t* hdr, bcf1_t* rec) {
        Out& o = out[i];
        o.unph = o.orphan = o.split = 0;
        bcf_unpack(rec, BCF_UN_ALL);
        // Only the phased record at a position is the one the GNN saw; other
        // records at the same position are written unchanged (issue #133).
        bool phased = false;
        {
            int32_t* gt=nullptr; int ng=0;
            if (bcf_get_genotypes(hdr,rec,&gt,&ng) > 0 && ng / nsmpl >= 2)
                phased = bcf_gt_is_phased(gt[si * (ng / nsmpl) + 1]);
            free(gt);
        }
        if (phased) {
            std::string ch = bcf_hdr_id2name(hdr, rec->rid);
            bool unph = false;
            auto pch = predictions_.find(ch);
            if (pch != predictions_.end()) {
                auto pit = pch->second.find(rec->pos);
                if (pit != pch->second.end() && shouldUnphase(pit->second)) unph = true;
            }
            if (!unph) {
                auto och = orphans_.find(ch);
                if (och != orphans_.end() && och->second.count(rec->pos)) {
                    unph = true;
                    o.orphan = 1;
                }
            }
            if (unph) {
                int32_t* gt=nullptr; int ng=0;
                if (bcf_get_genotypes(hdr,rec,&gt,&ng) > 0 && ng / nsmpl >= 2) {
                    // An unphased genotype carries no order, and the VCF spec
                    // writes it ascending, so sort the alleles rather than
                    // leaving "1/0" behind.
                    int32_t* sg = gt + si * (ng / nsmpl);
                    int a0 = bcf_gt_allele(sg[0]), a1 = bcf_gt_allele(sg[1]);
                    if (a0 > a1) std::swap(a0, a1);
                    sg[0]=bcf_gt_unphased(a0);
                    sg[1]=bcf_gt_unphased(a1);
                    bcf_update_genotypes(hdr,rec,gt,ng);
                }
                free(gt);
                setSamplePS(hdr, rec, si, bcf_int32_missing);
                o.unph = 1;
            } else if (params_.split_blocks) {
                // Variant stays phased but its block was split → new PS
                auto rch = ps_reassign_.find(ch);
                if (rch != ps_reassign_.end()) {
                    auto it = rch->second.find(rec->pos);
                    if (it != rch->second.end()) {
                        setSamplePS(hdr, rec, si, it->second);
                        o.split = 1;
                    }
                }
            }
        }
        o.text.l = 0;
        o.formatted = !(hdr->dirty && bcf_hdr_sync(hdr) < 0)
                      && bcf_hdr_nsamples(hdr) == rec->n_sample
                      && vcf_format(hdr, rec, &o.text) == 0;
    },
        [&](size_t i, bcf_hdr_t* hdr, bcf1_t* rec) {
        Out& o = out[i];
        nu += o.unph; norphan += o.orphan; nsplit += o.split;
        // bcf_write reports why the record could not be formatted
        if (!o.formatted && bcf_write(ofp, hdr, rec) == 0) return;
        if (!o.formatted || hwrite(ofp->fp.hfile, o.text.s, o.text.l) != (ssize_t)o.text.l) {
            std::cerr << "[GNN] ERROR: failed to write record to " << out_path << "\n";
            std::exit(EXIT_FAILURE);
        }
    });
    for (auto& o : out) free(o.text.s);
    bcf_hdr_destroy(hdr);
    hts_close(ifp);
    if (hts_close(ofp) != 0) {
        std::cerr << "[GNN] ERROR: failed to close " << out_path << "\n";
        std::exit(EXIT_FAILURE);
    }
    std::cerr << "[GNN] " << out_path << ": unphased " << nu
              << " (" << norphan << " left alone in their block)"
              << ", PS-split " << nsplit << " variants\n";
}

// ═══════════════════════════════════════════════════════════
void GNNModule::parseSecondaryVCF(const std::string& path, GnnVariantType vtype) {
    htsFile* fp = hts_open(path.c_str(), "r");
    // Fatal, like parseVCF: writeOutputVCF re-reads this file and would
    // fail anyway, only after the whole GNN pass.
    if (!fp) {
        std::cerr << "[GNN] ERROR: cannot open " << path << "\n";
        std::exit(EXIT_FAILURE);
    }
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        std::cerr << "[GNN] ERROR: cannot read VCF header of " << path << "\n";
        std::exit(EXIT_FAILURE);
    }
    int n = 0;
    // Run of consecutive heterozygous 5mC positions, as phase merges them
    // (a run is on one contig: same contig name, same header id)
    bool run_any = false;
    std::string run_chrom;
    int run_start = -1, run_last = -2;
    // each record of a batch: heterozygous, phased, and its variant
    std::vector<VariantInfo> got(kVcfBatch);
    std::vector<char> het(kVcfBatch), phased(kVcfBatch);
    readVcfRecords(fp, hdr, std::min(params_.threads, kVcfThreads), kVcfBatch,
        [&](size_t i, bcf_hdr_t* hdr, bcf1_t* rec) {
        het[i] = 0;
        bcf_unpack(rec, BCF_UN_ALL);
        int32_t* gt = nullptr; int ng = 0;
        if (bcf_get_genotypes(hdr, rec, &gt, &ng) < 0) return;
        if (ng < 2) { free(gt); return; }
        bool ph = bcf_gt_is_phased(gt[1]);
        int a0 = bcf_gt_allele(gt[0]), a1 = bcf_gt_allele(gt[1]);
        free(gt);
        if (a0 == a1) return;
        het[i] = 1;
        phased[i] = ph;
        VariantInfo& vi = got[i];
        vi = VariantInfo();
        vi.chrom = bcf_hdr_id2name(hdr, rec->rid);
        vi.pos = rec->pos;
        if (!ph) return;
        vi.gt_ref = a0; vi.gt_alt = a1;
        vi.is_phased = true; vi.pe = vi.h1 = vi.h2 = 0; vi.ps = -1;
        int rl = (int)strlen(rec->d.allele[0]);
        int al = rec->n_allele > 1 ? (int)strlen(rec->d.allele[1]) : rl;
        vi.indel_len = abs(al - rl);
        vi.vtype = vtype;
        float* fv = nullptr; int nv = 0;
        if (bcf_get_info_float(hdr, rec, "PE", &fv, &nv) > 0) { vi.pe = fv[0]; free(fv); fv=nullptr; nv=0; }
        if (bcf_get_info_float(hdr, rec, "H1", &fv, &nv) > 0) { vi.h1 = fv[0]; free(fv); fv=nullptr; nv=0; }
        if (bcf_get_info_float(hdr, rec, "H2", &fv, &nv) > 0) { vi.h2 = fv[0]; free(fv); fv=nullptr; nv=0; }
        int32_t* ps = nullptr; int np = 0;
        if (bcf_get_format_int32(hdr, rec, "PS", &ps, &np) > 0) { vi.ps = ps[0]; free(ps); }
        vi.hp_len = vi.gc_content = vi.str_context = vi.seq_entropy = 0;
    },
        [&](size_t i, bcf_hdr_t*, bcf1_t*) {
        if (!het[i]) return;
        VariantInfo& vi = got[i];
        int link_pos = -1;
        if (vtype == GVT_METHYL) {
            if (!run_any || vi.chrom != run_chrom || vi.pos != run_last + 1)
                run_start = vi.pos;
            if (!run_any || vi.chrom != run_chrom) run_chrom = vi.chrom;
            run_any = true;
            run_last = vi.pos;
            if (run_start != vi.pos) link_pos = run_start;
        }
        if (!phased[i]) return;
        vi.link_pos = link_pos;
        variants_[vi.chrom].push_back(vi);
        n++;
    });
    bcf_hdr_destroy(hdr); hts_close(fp);
    const char* tname[] = {"SNV","INDEL","SV","METHYL"};
    std::cerr << "  " << n << " phased " << tname[vtype] << " variants\n";
}

// ═══════════════════════════════════════════════════════════
// As in training: every symbol except N counts, IUPAC codes included.
float GNNModule::calcSeqEntropy(const std::string& s) {
    int c[256]={}, t=0; for(char ch:s){unsigned char u=std::toupper((unsigned char)ch); if(u!='N'){c[u]++;t++;}}
    if(!t) return 0;
    float e=0; for(int i=0;i<256;++i) if(c[i]>0){float p=(float)c[i]/t; e-=p*std::log2(p);} return e;
}
int GNNModule::calcHPLength(const std::string& s, int c) {
    if(c<0||c>=(int)s.size()) return 0;
    char b=std::toupper(s[c]); if(b=='N') return 0;
    int l=1; for(int i=c-1;i>=0&&std::toupper(s[i])==b;--i) l++;
    for(int i=c+1;i<(int)s.size()&&std::toupper(s[i])==b;++i) l++;
    return std::min(l,20);
}
// As in training (prepare_gnn_data.py): N counts in the denominator.
float GNNModule::calcGCContent(const std::string& s) {
    if (s.empty()) return 0.5f;
    int gc=0; for(char c:s){char u=std::toupper(c);if(u=='G'||u=='C')gc++;} return (float)gc/s.size();
}
float GNNModule::calcSTRContext(const std::string& s, int c) {
    if(s.size()<6||c<0||c>=(int)s.size()) return 0;
    int st=std::max(0,c-HP_RADIUS), en=std::min((int)s.size(),c+HP_RADIUS+1);
    const char* w=s.data()+st; const int wl=en-st; int best=0;
    // unit w[i, i+ul) without N, repeated back to back from i
    for(int ul=2;ul<=4;++ul) for(int i=0;i+ul<=wl;++i){
        if(std::memchr(w+i,'N',ul)) continue;
        int cnt=0; for(int j=i;j+ul<=wl;j+=ul) if(std::memcmp(w+j,w+i,ul)==0)cnt++;else break;
        if(cnt>=2&&cnt>best)best=cnt;}
    return std::min((float)best/10.0f,1.0f);
}
