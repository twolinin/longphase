/*
 * GNNProcess.h - LongPhase GNN post-hoc phasing correction module
 */
#ifndef GNN_PROCESS_H
#define GNN_PROCESS_H

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <cmath>
#include <mutex>
#include <onnxruntime_cxx_api.h>

constexpr int NODE_FEAT_DIM = 31;
constexpr int EDGE_FEAT_DIM = 7;
constexpr int HP_RADIUS     = 15;
constexpr int GC_RADIUS     = 50;
constexpr int NEIGH_RADIUS  = 500;
constexpr int MAX_NODES     = 256;

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

// DOT edge: 1-indexed alleles, 0-based positions
struct DotEdge {
    int src_pos, src_allele, dst_pos, dst_allele;
    float weight;
};

struct Prediction {
    float prob_error = 0.0f;
    bool  is_bridge  = false;
    int   count      = 0;
};

class GNNModule {
public:
    struct Params {
        std::string model_path, vcf_path, dot_prefix, reference_path, output_vcf;
        std::string sv_vcf, mod_vcf, output_sv_vcf, output_mod_vcf;
        float break_threshold = 0.30f;
        float pe_threshold    = 0.80f;
        int   window          = 20;
        int   threads         = 4;
        bool  respect_bridge  = false;
    };

    explicit GNNModule(const Params& p);
    ~GNNModule();
    int run();

private:
    Params params_;
    Ort::Env                      ort_env_;
    std::unique_ptr<Ort::Session> ort_session_;
    Ort::SessionOptions           ort_options_;
    Ort::MemoryInfo               memory_info_;

    std::map<std::string, std::vector<VariantInfo>> variants_;
    // dot_edges_[chrom][src_pos] = list of edges from that position
    std::map<std::string, std::unordered_map<int, std::vector<DotEdge>>> dot_edges_;

    std::mutex pred_mutex_;
    std::map<std::string, std::unordered_map<int, Prediction>> predictions_;

    void loadModel();
    void parseVCF();
    void parseSecondaryVCF(const std::string& path, GnnVariantType vtype);
    void parseDotFiles();
    void computeGenomicFeatures();
    // Pre-computed per-chromosome data
    struct ChromData {
        std::unordered_map<int, int> ps_count;
        std::vector<int> all_pos;
        std::vector<int> indel_pos;
        std::unordered_set<int> bridge_set;   // precomputed bridge vertices
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

#endif
