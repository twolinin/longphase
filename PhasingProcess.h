#ifndef PHASINGPROCESS_H
#define PHASINGPROCESS_H
#include "Util.h"

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include <atomic>
#include <memory>



struct PhasingParameters
{
    int numThreads;
    int distance;
    std::string snpFile;
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
struct DotEdge {
    int src_pos, src_allele, dst_pos, dst_allele;
    float weight;
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
        float break_threshold = 0.30f;
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
    void setDotEdges(std::map<std::string, std::vector<DotEdge>>&& edges);
    int run();

private:
    Params params_;

    // Holds the decoded weights. forward() keeps no mutable state, so a
    // single const instance is shared by every worker thread.
    std::unique_ptr<const gnn::Model> model_;

    std::map<std::string, std::vector<VariantInfo>> variants_;
    // dot_edges_[chrom][src_pos] = list of edges from that position
    std::map<std::string, std::unordered_map<int, std::vector<DotEdge>>> dot_edges_;
    std::map<std::string, std::vector<DotEdge>> pending_edges_;

    std::mutex pred_mutex_;
    // Windows skipped because they exceed MAX_NODES
    std::atomic<int> skipped_windows_{0};
    std::map<std::string, std::unordered_map<int, Prediction>> predictions_;

    // Block splits: [chrom][pos] = new_ps_id (for variants whose block
    // was split after a bridge was unphased)
    std::map<std::string, std::unordered_map<int, int>> ps_reassign_;

    void loadModel();
    void parseVCF();
    void parseSecondaryVCF(const std::string& path, GnnVariantType vtype);
    void loadDotEdges();
    void computeGenomicFeatures();
    void computeBlockSplits();
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
