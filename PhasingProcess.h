#ifndef PHASINGPROCESS_H
#define PHASINGPROCESS_H
#include "Util.h"

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include <atomic>
#include <memory>
#include <htslib/vcf.h>



struct PhasingParameters
{
    int numThreads;
    int distance;
    std::string snpFile;
    // --sample: the sample of a multi-sample VCF to phase (empty: the first)
    std::string sample;
    std::string svFile;
    std::vector<std::string> bamFile;
    std::string modFile="";
    std::string fastaFile;
    std::string resultPrefix;
    bool generateDot;
    bool isONT;
    bool isPB;
    bool phaseIndel;
    int indelQuality;  // Indel quality filter threshold, default is 0 (disabled)
    
    int connectAdjacent;
    int mappingQuality;
    double mismatchRate;
    
    int baseQuality;
    double edgeWeight;
    
    double snpConfidence;
    double readConfidence;
    
    double edgeThreshold;
    double overlapThreshold;

    int svWindow;
    double svThreshold; 

    // GNN stage of phasing (off with --disableGNN)
    bool enableGNN;
    float gnnBreakThreshold;
    float gnnPeThreshold;
    int gnnWindow;
    bool gnnRespectBridge;
    bool gnnSplitBlocks;
    
    std::string version;
    std::string command;
};

// ─── GNN ──────────────────────────────────────────────────
// The GNN stage of phasing (off with --disableGNN). The model and its
// compiled-in weights are only included by PhasingProcess.cpp.
namespace gnn { class Model; }

// One edge of the phasing graph: 1-indexed alleles, 0-based positions.
// The source position is the one the edge is listed under in GnnEdgeList.
struct GnnEdge {
    int dst_pos;
    float weight;
    int8_t src_allele, dst_allele;
};

// The graph edges of one chromosome, grouped by source position. The edges
// of src[i] are edges[start[i]] up to the start of the next group, in the
// order they were added.
struct GnnEdgeList {
    std::vector<int> src;
    std::vector<size_t> start;
    std::vector<GnnEdge> edges;

    struct Range {
        const GnnEdge *b, *e;
        const GnnEdge *begin() const { return b; }
        const GnnEdge *end() const { return e; }
    };

    void add(int srcPos, int srcAllele, int dstPos, int dstAllele, float weight) {
        if (src.empty() || src.back() != srcPos) {
            if (!src.empty() && srcPos < src.back())
                sorted = false;
            src.push_back(srcPos);
            start.push_back(edges.size());
        }
        edges.push_back({dstPos, weight, (int8_t)srcAllele, (int8_t)dstAllele});
    }
    void swap(GnnEdgeList &o) {
        src.swap(o.src);
        start.swap(o.start);
        edges.swap(o.edges);
        std::swap(sorted, o.sorted);
    }
    // Group the edges by ascending source position. phase adds them in that
    // order already; otherwise sort them, keeping the order within a source.
    void finish();
    size_t groups() const { return src.size(); }
    Range group(size_t i) const {
        size_t e = i + 1 < start.size() ? start[i + 1] : edges.size();
        return {edges.data() + start[i], edges.data() + e};
    }
    // 1 if srcPos has edges
    size_t count(int srcPos) const {
        auto it = std::lower_bound(src.begin(), src.end(), srcPos);
        return it != src.end() && *it == srcPos;
    }
    // the edges of srcPos, which must have some
    Range at(int srcPos) const {
        return group(std::lower_bound(src.begin(), src.end(), srcPos) - src.begin());
    }

private:
    bool sorted = true;
};

enum GnnVariantType { GVT_SNV=0, GVT_INDEL, GVT_SV, GVT_METHYL };

struct VariantInfo {
    std::string chrom;
    int   pos;           // 0-based
    float pe;
    int   gt_ref, gt_alt; // genotype allele digits
    float h1, h2;         // VCF INFO H1/H2
    int   ps;
    bool  is_phased;
    int   indel_len;      // |len(alt) - len(ref)|
    float hp_len, gc_content, str_context, seq_entropy;
    GnnVariantType vtype = GVT_SNV;
    // 5mC only: phase merges a run of consecutive 5mC positions (the two
    // strands of a CpG) into one graph node at the first position of the
    // run. Position of that node when it is another site, otherwise -1.
    int   link_pos = -1;
};

struct Prediction {
    float prob_error = 0.0f;
    bool  is_bridge  = false;
    int   count      = 0;
};

class GNNModule {
public:
    struct Params {
        std::string vcf_path, reference_path, output_vcf;
        std::string sv_vcf, mod_vcf, output_sv_vcf, output_mod_vcf;
        // --sample; a VCF without this sample (SV, 5mC) uses its first sample
        std::string sample;
        float break_threshold = 0.35f;
        float pe_threshold    = 0.80f;
        int   window          = 20;
        int   threads         = 4;
        bool  respect_bridge  = false;
        bool  split_blocks    = true;
    };

    explicit GNNModule(const Params& p);
    ~GNNModule();
    // Graph edges from phase, per chromosome, in the order phase emitted
    // them. Taken over by run().
    void setDotEdges(std::map<std::string, GnnEdgeList>&& edges);
    int run();

private:
    Params params_;

    // Holds the decoded weights. forward() keeps no mutable state, so a
    // single const instance is shared by every worker thread.
    std::unique_ptr<const gnn::Model> model_;

    std::map<std::string, std::vector<VariantInfo>> variants_;
    // dot_edges_[chrom]: the edges of each source position
    std::map<std::string, GnnEdgeList> dot_edges_;
    std::map<std::string, GnnEdgeList> pending_edges_;

    std::mutex pred_mutex_;
    // Windows skipped because they exceed MAX_NODES
    std::atomic<int> skipped_windows_{0};
    std::map<std::string, std::unordered_map<int, Prediction>> predictions_;

    // Block splits: [chrom][pos] = new_ps_id (for variants whose block
    // was split after a bridge was unphased)
    std::map<std::string, std::unordered_map<int, int>> ps_reassign_;

    // Variants unphased only because GNN correction left them alone in
    // their phase set: [chrom] = positions
    std::map<std::string, std::unordered_set<int>> orphans_;

    void loadModel();
    void parseVCF();
    void parseSecondaryVCF(const std::string& path, GnnVariantType vtype);
    void loadDotEdges();
    void computeGenomicFeatures();
    void computeBlockSplits();
    void computeOrphans();
    void shareMergedPredictions();
    bool shouldUnphase(const Prediction& pr) const;
    // Pre-computed per-chromosome data
    struct ChromData {
        std::unordered_map<int, int> ps_count;
        std::vector<int> all_pos;
        std::vector<int> indel_pos;
        std::unordered_map<int, int> pos_to_idx; // pos → index in variant list
    };

    void processChromosome(const std::string& chrom, int trigger_idx,
                           const ChromData& cd);
    void writeOutputVCF();
    void writeOutputVCF(const std::string& in_path, const std::string& out_path);
    int  sampleIndex(bcf_hdr_t* hdr) const;
    static void setSamplePS(bcf_hdr_t* hdr, bcf1_t* rec, int si, int32_t ps);

    // Feature computation helpers (static, matching prepare_gnn_data.py)
    static float calcSeqEntropy(const std::string& seq);
    static int   calcHPLength(const std::string& seq, int center);
    static float calcGCContent(const std::string& seq);
    static float calcSTRContext(const std::string& seq, int center);
};

class PhasingProcess
{

    public:
        PhasingProcess(PhasingParameters params);
        ~PhasingProcess();

};


#endif
