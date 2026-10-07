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
#include <thread>
#include <set>
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

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
    // Graph edges for the GNN, per chromosome; pre-created for the
    // same reason as chrPhasingResult.
    std::map<std::string, std::vector<DotEdge> > chrGnnEdges;
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
        std::string &chr_reference = fastaParser.chrString.at(*chrIter);
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
        clip->getCNVInterval(clipCount);
        

        // create a graph object and prepare to phasing.
        VairiantGraph *vGraph = new VairiantGraph(chr_reference, params, (*chrIter));
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
        std::map<std::string, std::string>().swap(fastaParser.chrString);
        PhasingResult().swap(mergedPhasingResult);
        ChrPhasingResult().swap(chrPhasingResult);

        begin = time(NULL);
        std::cerr<< "\nGNN ...\n";
        GNNModule::Params gnnParams;
        gnnParams.vcf_path        = phasedPrefix + ".vcf";
        gnnParams.reference_path  = params.fastaFile;
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
    bcf1_t* rec = bcf_init();
    int n = 0;
    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);
        int32_t* gt = nullptr; int ng = 0;
        if (bcf_get_genotypes(hdr, rec, &gt, &ng) < 0) continue;
        if (ng < 2) { free(gt); continue; }
        if (!bcf_gt_is_phased(gt[1])) { free(gt); continue; }
        int a0 = bcf_gt_allele(gt[0]), a1 = bcf_gt_allele(gt[1]);
        free(gt);
        if (a0 == a1) continue;

        VariantInfo vi;
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
        if (bcf_get_format_int32(hdr, rec, "PS", &ps, &np) > 0) { vi.ps = ps[0]; free(ps); }
        vi.hp_len = vi.gc_content = vi.str_context = vi.seq_entropy = 0;
        vi.vtype = (vi.indel_len == 0) ? GVT_SNV : GVT_INDEL;
        variants_[vi.chrom].push_back(vi);
        n++;
    }
    bcf_destroy(rec); bcf_hdr_destroy(hdr); hts_close(fp);
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
void GNNModule::setDotEdges(std::map<std::string, std::vector<DotEdge>>&& edges) {
    pending_edges_ = std::move(edges);
}

void GNNModule::loadDotEdges() {
    // Collect chromosome list and pre-initialize dot_edges_ map
    std::vector<std::string> chroms;
    for (auto& kv : variants_) {
        const std::string& ch = kv.first;
        chroms.push_back(ch);
        dot_edges_[ch]; // pre-create entry to avoid concurrent map insertion
    }

    // Index all chromosomes in parallel
    std::atomic<int> total{0};
    std::vector<std::thread> threads;
    std::atomic<size_t> cidx{0};

    int n_threads = std::min(params_.threads, (int)chroms.size());
    for (int t = 0; t < n_threads; ++t) {
        threads.emplace_back([&]() {
            while (true) {
                size_t i = cidx.fetch_add(1);
                if (i >= chroms.size()) return;
                auto& ch = chroms[i];
                auto pit = pending_edges_.find(ch);
                if (pit == pending_edges_.end()) continue;
                const std::vector<DotEdge>& edges = pit->second;

                // phase emits all edges of a source position together, so
                // store each run at its exact size.
                auto& em = dot_edges_[ch];
                size_t r = 0;
                while (r < edges.size()) {
                    size_t e = r + 1;
                    while (e < edges.size() && edges[e].src_pos == edges[r].src_pos) ++e;
                    auto& v = em[edges[r].src_pos];
                    if (v.empty()) v.assign(edges.begin() + r, edges.begin() + e);
                    else v.insert(v.end(), edges.begin() + r, edges.begin() + e);
                    r = e;
                }
                total += (int)edges.size();
                std::vector<DotEdge>().swap(pit->second);
            }
        });
    }
    for (auto& t : threads) t.join();
    pending_edges_.clear();
    std::cerr << "  " << total.load() << " graph edges\n";
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
            // Each thread gets its own faidx handle
            faidx_t* fai = fai_load(params_.reference_path.c_str());
            if (!fai) return;
            while (true) {
                size_t i = cidx.fetch_add(1);
                if (i >= chroms.size()) { fai_destroy(fai); return; }
                auto& ch = chroms[i];
                for (auto& v : variants_[ch]) {
                    // The model was trained with these four features set
                    // to 0 for SVs, whose breakpoint context is not
                    // meaningful, so leave them at 0 here too.
                    if (v.vtype == GVT_SV) continue;
                    int len = 0;
                    int s1 = std::max(0, v.pos - HP_RADIUS);
                    char* seq = faidx_fetch_seq(fai, ch.c_str(), s1, v.pos+HP_RADIUS, &len);
                    if (seq && len > 0) {
                        std::string ss(seq, len);
                        v.seq_entropy = calcSeqEntropy(ss);
                        v.hp_len = (float)calcHPLength(ss, v.pos-s1);
                        v.str_context = calcSTRContext(ss, v.pos-s1);
                        free(seq);
                    }
                    int s2 = std::max(0, v.pos - GC_RADIUS);
                    seq = faidx_fetch_seq(fai, ch.c_str(), s2, v.pos+GC_RADIUS, &len);
                    if (seq && len > 0) {
                        v.gc_content = calcGCContent(std::string(seq, len));
                        free(seq);
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
    std::unordered_set<int> wpos_set;
    std::unordered_map<int, const VariantInfo*> pos_to_var;
    for (int i = wi_start; i <= wi_end; ++i) {
        wvars.push_back(&vlist[i]);
        wpos_set.insert(vlist[i].pos);
        pos_to_var[vlist[i].pos] = &vlist[i];
    }
    int n_var = (int)wvars.size();
    int N = n_var * 2;
    if (N > MAX_NODES) { skipped_windows_.fetch_add(1); return; }

    float max_offset = 1.0f;
    for (auto* v : wvars) {
        float d = std::abs((float)(v->pos - center_pos));
        if (d > max_offset) max_offset = d;
    }

    // h1/h2 from DOT first edge per allele
    std::unordered_map<int, float> dot_h1, dot_h2;
    for (auto* v : wvars) {
        if (!edge_map.count(v->pos)) continue;
        for (auto& e : edge_map.at(v->pos)) {
            if (e.src_allele == 1 && !dot_h1.count(v->pos)) dot_h1[v->pos] = e.weight;
            else if (e.src_allele == 2 && !dot_h2.count(v->pos)) dot_h2[v->pos] = e.weight;
        }
    }

    // out_degree, out_weight_std
    std::unordered_map<int, float> out_deg, out_wstd;
    for (auto* v : wvars) {
        if (!edge_map.count(v->pos)) { out_deg[v->pos]=0; out_wstd[v->pos]=0; continue; }
        std::set<int> dsts; std::vector<float> ws;
        for (auto& e : edge_map.at(v->pos)) { dsts.insert(e.dst_pos); ws.push_back(e.weight); }
        out_deg[v->pos] = std::min((float)dsts.size()/10.0f, 1.0f);
        if (ws.size() > 1) {
            float mu=0; for(float w:ws) mu+=w; mu/=ws.size();
            float var=0; for(float w:ws) var+=(w-mu)*(w-mu); var/=ws.size();
            out_wstd[v->pos] = std::min(std::sqrt(var)/100.0f, 1.0f);
        } else out_wstd[v->pos] = 0;
    }

    // window_mean_votes (VCF h1+h2, full chromosome ±10)
    std::unordered_map<int, float> wmean;
    for (int i = wi_start; i <= wi_end; ++i) {
        int lo = std::max(0, i-10), hi = std::min((int)vlist.size()-1, i+10);
        float sum=0; int cnt=0;
        for (int j=lo; j<=hi; ++j) { sum += vlist[j].h1 + vlist[j].h2; cnt++; }
        wmean[vlist[i].pos] = cnt > 0 ? sum/cnt : 1.0f;
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
    std::unordered_map<int, std::unordered_set<int>> pos_adj;
    for (int pos : wpos_set) {
        if (!edge_map.count(pos)) continue;
        for (auto& e : edge_map.at(pos)) {
            if (!wpos_set.count(e.dst_pos)) continue;
            if (e.src_pos == e.dst_pos) continue;
            pos_adj[e.src_pos].insert(e.dst_pos);
            pos_adj[e.dst_pos].insert(e.src_pos);
        }
    }

    // Connected components of wpos_set, optionally with one node removed.
    auto n_components = [&](int exclude) -> int {
        std::unordered_set<int> vis;
        int count = 0;
        for (int start : wpos_set) {
            if (start == exclude) continue;
            if (vis.count(start)) continue;
            ++count;
            std::vector<int> q{start};
            while (!q.empty()) {
                int node = q.back(); q.pop_back();
                if (vis.count(node)) continue;
                vis.insert(node);
                auto it = pos_adj.find(node);
                if (it == pos_adj.end()) continue;
                for (int nb : it->second)
                    if (nb != exclude && wpos_set.count(nb) && !vis.count(nb))
                        q.push_back(nb);
            }
        }
        return count;
    };

    const int comp_before = n_components(-1);
    std::unordered_set<int> bridge_set;
    for (int pos : wpos_set) {
        if (wpos_set.size() <= 1) break;
        auto it = pos_adj.find(pos);
        if (it == pos_adj.end()) continue;
        // Needs at least one neighbour that survives the removal.
        bool has_active = false;
        for (int nb : it->second)
            if (nb != pos && wpos_set.count(nb)) { has_active = true; break; }
        if (!has_active) continue;
        if (n_components(pos) > comp_before) bridge_set.insert(pos);
    }

    // out_wsum for weight_ratio
    std::unordered_map<int64_t, float> out_wsum;
    auto make_key = [](int pos, int al) -> int64_t { return ((int64_t)pos << 4) | al; };
    for (int pos : wpos_set) {
        if (!edge_map.count(pos)) continue;
        for (auto& e : edge_map.at(pos)) {
            if (!wpos_set.count(e.dst_pos)) continue;
            out_wsum[make_key(e.src_pos, e.src_allele)] += e.weight;
        }
    }

    // ── Node features [N × NODE_FEAT_DIM] ──
    std::vector<float> nf(N * NODE_FEAT_DIM, 0.0f);
    std::unordered_map<int64_t, int> pa_to_nid;

    for (int i = 0; i < n_var; ++i) {
        const auto* v = wvars[i];
        float h1v = dot_h1.count(v->pos) ? dot_h1[v->pos] : v->h1;
        float h2v = dot_h2.count(v->pos) ? dot_h2[v->pos] : v->h2;
        float total = h1v + h2v;
        float vote_r = total > 0 ? std::max(h1v,h2v)/total : 0.0f;
        float log_tot = std::log2(1.0f+total);
        float lcd = std::log10(1.0f + std::abs((float)(v->pos-center_pos)));
        float ps_blk = std::log1p((float)(ps_count.count(v->ps)?ps_count.at(v->ps):1))
                       / std::log1p(100.0f);
        float wm = wmean.count(v->pos) ? wmean[v->pos] : 1.0f;
        float rvd = std::log2(total / (wm+EPS) + EPS);
        float bridge_f = bridge_set.count(v->pos) ? 1.0f : 0.0f;

        for (int a = 0; a < 2; ++a) {
            int allele = a+1, nid = i*2+a;
            pa_to_nid[make_key(v->pos, allele)] = nid;
            int hp_val = (allele==1) ? v->gt_ref : v->gt_alt;
            float* f = &nf[nid * NODE_FEAT_DIM];
            f[0]=v->pe; f[1]=(float)hp_val; f[2]=(float)v->gt_ref; f[3]=(float)v->gt_alt;
            f[4]=(v->pos==center_pos)?1.0f:0.0f;
            f[5]=(v->ps==center_ps && v->ps!=-1)?1.0f:0.0f;
            f[6]=(float)a; f[7]=(float)(v->pos-center_pos)/max_offset;
            f[8]=h1v; f[9]=h2v; f[10]=total; f[11]=vote_r;
            f[12]=log_tot; f[13]=std::log1p(std::max(0.0f,h1v)); f[14]=std::log1p(std::max(0.0f,h2v));
            f[15]=v->hp_len; f[16]=v->gc_content; f[17]=v->str_context; f[18]=v->seq_entropy;
            f[19]=neigh_dens(v->pos); f[20]=out_deg[v->pos]; f[21]=ps_blk;
            f[22]=lcd; f[23]=dist_indel(v->pos); f[24]=rvd;
            f[25]=out_wstd[v->pos]; f[26]=bridge_f;
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

    for (int pos : wpos_set) {
        if (!edge_map.count(pos)) continue;
        for (auto& e : edge_map.at(pos)) {
            if (!wpos_set.count(e.dst_pos)) continue;
            auto sit = pa_to_nid.find(make_key(e.src_pos, e.src_allele));
            auto dit = pa_to_nid.find(make_key(e.dst_pos, e.dst_allele));
            if (sit == pa_to_nid.end() || dit == pa_to_nid.end()) continue;
            int sn = sit->second, dn = dit->second;
            float same_hp = (e.src_allele == e.dst_allele) ? 1.0f : 0.0f;
            auto* sv = pos_to_var[e.src_pos]; auto* dv = pos_to_var[e.dst_pos];
            float inter_ps = (sv && dv && sv->ps != dv->ps) ? 1.0f : 0.0f;
            float wsum = out_wsum[make_key(e.src_pos, e.src_allele)];
            int row = dn, col = sn;
            adj[row*N+col] = 1.0f;
            float* ep = &ef[(row*N+col)*EDGE_FEAT_DIM];
            ep[0] = 1.0f/(1.0f+std::exp(-e.weight));
            ep[1] = std::log10(1.0f+std::abs((float)(e.dst_pos-e.src_pos)));
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
        bool br = bridge_set.count(wvars[i]->pos) > 0;
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

            // Build adjacency (only same-PS, both-phased edges)
            std::unordered_map<int, std::vector<int>> adj;
            for (auto& kv_src : em) {
                const auto& src = kv_src.first;
                auto& elist = kv_src.second;
                if (unphase_pos.count(src)) continue;
                auto it_s = pos_to_ps.find(src);
                if (it_s == pos_to_ps.end()) continue;
                for (auto& e : elist) {
                    int dst = e.dst_pos;
                    if (unphase_pos.count(dst)) continue;
                    auto it_d = pos_to_ps.find(dst);
                    if (it_d == pos_to_ps.end()) continue;
                    if (it_s->second != it_d->second) continue;
                    adj[src].push_back(dst);
                    adj[dst].push_back(src);
                }
            }
            // A 5mC site merged into another site's graph node has no edges
            // of its own; it is connected through that node.
            for (auto& v : vl) {
                if (v.link_pos < 0) continue;
                auto it_s = pos_to_ps.find(v.pos);
                auto it_d = pos_to_ps.find(v.link_pos);
                if (it_s == pos_to_ps.end() || it_d == pos_to_ps.end()) continue;
                if (it_s->second != it_d->second) continue;
                adj[v.pos].push_back(v.link_pos);
                adj[v.link_pos].push_back(v.pos);
            }

            std::unordered_map<int, int> local_reassign;
            int local_splits = 0, local_new = 0;

            // Per-PS connected components via BFS
            for (auto& kv_ps : ps_blocks) {
                auto& positions = kv_ps.second;
                if (positions.size() <= 1) continue;

                std::unordered_set<int> pos_set(positions.begin(), positions.end());
                std::unordered_set<int> visited;
                std::vector<std::vector<int>> components;

                for (int start : positions) {
                    if (visited.count(start)) continue;
                    std::vector<int> comp;
                    std::vector<int> queue = {start};
                    while (!queue.empty()) {
                        int node = queue.back(); queue.pop_back();
                        if (visited.count(node)) continue;
                        visited.insert(node);
                        comp.push_back(node);
                        auto it = adj.find(node);
                        if (it != adj.end())
                            for (int nb : it->second)
                                if (!visited.count(nb) && pos_set.count(nb))
                                    queue.push_back(nb);
                    }
                    components.push_back(std::move(comp));
                }

                if (components.size() <= 1) continue;

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
    bcf1_t* rec = bcf_init();
    int nu = 0, nsplit = 0, norphan = 0;
    while (bcf_read(ifp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);
        std::string ch = bcf_hdr_id2name(hdr, rec->rid);
        bool unph = false;
        if (predictions_.count(ch) && predictions_[ch].count(rec->pos)) {
            auto& pr = predictions_[ch][rec->pos];
            if (shouldUnphase(pr)) unph = true;
        }
        if (!unph && orphans_.count(ch) && orphans_[ch].count(rec->pos)) {
            unph = true;
            norphan++;
        }
        if (unph) {
            int32_t* gt=nullptr; int ng=0;
            if (bcf_get_genotypes(hdr,rec,&gt,&ng)>=2) {
                // An unphased genotype carries no order, and the VCF spec
                // writes it ascending, so sort the alleles rather than
                // leaving "1/0" behind.
                int a0 = bcf_gt_allele(gt[0]), a1 = bcf_gt_allele(gt[1]);
                if (a0 > a1) std::swap(a0, a1);
                gt[0]=bcf_gt_unphased(a0);
                gt[1]=bcf_gt_unphased(a1);
                bcf_update_genotypes(hdr,rec,gt,ng);
            }
            free(gt);
            bcf_update_format_int32(hdr,rec,"PS",nullptr,0);
            nu++;
        } else if (params_.split_blocks && ps_reassign_.count(ch)) {
            // Variant stays phased but its block was split → new PS
            auto& rmap = ps_reassign_[ch];
            auto it = rmap.find(rec->pos);
            if (it != rmap.end()) {
                int32_t new_ps = it->second;
                bcf_update_format_int32(hdr, rec, "PS", &new_ps, 1);
                nsplit++;
            }
        }
        if (bcf_write(ofp, hdr, rec) < 0) {
            std::cerr << "[GNN] ERROR: failed to write record to " << out_path << "\n";
            std::exit(EXIT_FAILURE);
        }
    }
    bcf_destroy(rec); bcf_hdr_destroy(hdr);
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
    bcf1_t* rec = bcf_init();
    int n = 0;
    // Run of consecutive heterozygous 5mC positions, as phase merges them
    int run_rid = -1, run_start = -1, run_last = -2;
    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);
        int32_t* gt = nullptr; int ng = 0;
        if (bcf_get_genotypes(hdr, rec, &gt, &ng) < 0) continue;
        if (ng < 2) { free(gt); continue; }
        bool phased = bcf_gt_is_phased(gt[1]);
        int a0 = bcf_gt_allele(gt[0]), a1 = bcf_gt_allele(gt[1]);
        free(gt);
        if (a0 == a1) continue;
        int link_pos = -1;
        if (vtype == GVT_METHYL) {
            if (rec->rid != run_rid || rec->pos != run_last + 1)
                run_start = (int)rec->pos;
            run_rid = rec->rid;
            run_last = (int)rec->pos;
            if (run_start != rec->pos) link_pos = run_start;
        }
        if (!phased) continue;
        VariantInfo vi;
        vi.chrom = bcf_hdr_id2name(hdr, rec->rid);
        vi.pos = rec->pos; vi.gt_ref = a0; vi.gt_alt = a1;
        vi.is_phased = true; vi.pe = vi.h1 = vi.h2 = 0; vi.ps = -1;
        int rl = (int)strlen(rec->d.allele[0]);
        int al = rec->n_allele > 1 ? (int)strlen(rec->d.allele[1]) : rl;
        vi.indel_len = abs(al - rl);
        vi.vtype = vtype;
        vi.link_pos = link_pos;
        float* fv = nullptr; int nv = 0;
        if (bcf_get_info_float(hdr, rec, "PE", &fv, &nv) > 0) { vi.pe = fv[0]; free(fv); fv=nullptr; nv=0; }
        if (bcf_get_info_float(hdr, rec, "H1", &fv, &nv) > 0) { vi.h1 = fv[0]; free(fv); fv=nullptr; nv=0; }
        if (bcf_get_info_float(hdr, rec, "H2", &fv, &nv) > 0) { vi.h2 = fv[0]; free(fv); fv=nullptr; nv=0; }
        int32_t* ps = nullptr; int np = 0;
        if (bcf_get_format_int32(hdr, rec, "PS", &ps, &np) > 0) { vi.ps = ps[0]; free(ps); }
        vi.hp_len = vi.gc_content = vi.str_context = vi.seq_entropy = 0;
        variants_[vi.chrom].push_back(vi);
        n++;
    }
    bcf_destroy(rec); bcf_hdr_destroy(hdr); hts_close(fp);
    const char* tname[] = {"SNV","INDEL","SV","METHYL"};
    std::cerr << "  " << n << " phased " << tname[vtype] << " variants\n";
}

// ═══════════════════════════════════════════════════════════
float GNNModule::calcSeqEntropy(const std::string& s) {
    int c[4]={}; for(char ch:s) switch(std::toupper(ch)){case'A':c[0]++;break;case'C':c[1]++;break;case'G':c[2]++;break;case'T':c[3]++;break;}
    int t=c[0]+c[1]+c[2]+c[3]; if(!t) return 0;
    float e=0; for(int i=0;i<4;++i) if(c[i]>0){float p=(float)c[i]/t; e-=p*std::log2(p);} return e;
}
int GNNModule::calcHPLength(const std::string& s, int c) {
    if(c<0||c>=(int)s.size()) return 0;
    char b=std::toupper(s[c]); if(b=='N') return 0;
    int l=1; for(int i=c-1;i>=0&&std::toupper(s[i])==b;--i) l++;
    for(int i=c+1;i<(int)s.size()&&std::toupper(s[i])==b;++i) l++;
    return std::min(l,20);
}
float GNNModule::calcGCContent(const std::string& s) {
    int gc=0,t=0; for(char c:s){char u=std::toupper(c);if(u=='G'||u=='C')gc++;if(u!='N')t++;} return t>0?(float)gc/t:0.5f;
}
float GNNModule::calcSTRContext(const std::string& s, int c) {
    if(s.size()<6||c<0||c>=(int)s.size()) return 0;
    int st=std::max(0,c-HP_RADIUS), en=std::min((int)s.size(),c+HP_RADIUS+1);
    std::string w=s.substr(st,en-st); int best=0;
    for(int ul=2;ul<=4;++ul) for(int i=0;i+ul<=(int)w.size();++i){
        std::string u=w.substr(i,ul); if(u.find('N')!=std::string::npos) continue;
        int cnt=0; for(int j=i;j+ul<=(int)w.size();j+=ul) if(w.substr(j,ul)==u)cnt++;else break;
        if(cnt>=2&&cnt>best)best=cnt;}
    return std::min((float)best/10.0f,1.0f);
}
