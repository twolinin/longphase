/*
 * gnn_module.cpp - LongPhase GNN correction (optimized)
 *
 * Optimizations vs initial version:
 *   1. DOT parsing: sscanf (was std::regex, 10-50x faster)
 *   2. Bridge detection: precomputed per chromosome (was per-window BFS)
 *   3. Variant lookup: unordered_map (was linear scan)
 *   4. Per-chromosome data precomputed once
 *   5. Removed debug output
 */
#include "GNNProcess.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <cstdio>
#include <thread>
#include <atomic>
#include <set>
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

static constexpr float EPS = 1e-6f;

GNNModule::GNNModule(const Params& p)
    : params_(p),
      ort_env_(ORT_LOGGING_LEVEL_WARNING, "longphase_gnn"),
      memory_info_(Ort::MemoryInfo::CreateCpu(
          OrtAllocatorType::OrtArenaAllocator, OrtMemTypeDefault))
{}
GNNModule::~GNNModule() = default;

// ═══════════════════════════════════════════════════════════
int GNNModule::run() {
    std::cerr << "[GNN] Loading ONNX model\n";
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
        for (auto& [ch, vl] : variants_)
            std::sort(vl.begin(), vl.end(),
                      [](const VariantInfo& a, const VariantInfo& b){ return a.pos < b.pos; });
    }
    std::cerr << "[GNN] Parsing DOT files\n";
    parseDotFiles();
    if (!params_.reference_path.empty()) {
        std::cerr << "[GNN] Computing genomic features\n";
        computeGenomicFeatures();
    }

    // Pre-compute per-chromosome data
    std::map<std::string, ChromData> chrom_data;
    std::vector<std::string> chroms;
    for (auto& [ch, vl] : variants_) chroms.push_back(ch);

    for (auto& ch : chroms) {
        auto& cd = chrom_data[ch];
        auto& vl = variants_[ch];
        for (auto& v : vl) cd.ps_count[v.ps]++;
        cd.all_pos.reserve(vl.size());
        for (auto& v : vl) cd.all_pos.push_back(v.pos);
        for (auto& v : vl) if (v.indel_len > 0) cd.indel_pos.push_back(v.pos);
        std::sort(cd.indel_pos.begin(), cd.indel_pos.end());

        // Precompute bridge vertices using full chromosome adjacency
        auto& em = dot_edges_[ch];
        std::unordered_map<int, std::unordered_set<int>> pos_adj;
        std::unordered_set<int> all_ps_set;
        for (auto& v : vl) {
            all_ps_set.insert(v.pos);
            if (!em.count(v.pos)) continue;
            for (auto& e : em.at(v.pos)) {
                pos_adj[v.pos].insert(e.dst_pos);
                pos_adj[e.dst_pos].insert(v.pos);
            }
        }
        // BFS for bridge detection on the full graph
        // A variant is a bridge if removing it increases connected components
        // Only check variants with PE >= threshold (potential triggers and their neighbors)
        for (auto& v : vl) {
            if (v.pe < params_.pe_threshold * 0.5f) continue; // only check relevant variants
            if (pos_adj.find(v.pos) == pos_adj.end()) continue;
            // Quick check: if degree <= 1, not a bridge
            if (pos_adj[v.pos].size() <= 1) continue;

            // Count components in local neighborhood with and without
            auto local = pos_adj[v.pos]; // neighbors
            local.insert(v.pos);
            // Add 2nd-hop neighbors
            std::unordered_set<int> extended;
            for (int p : local) {
                extended.insert(p);
                if (pos_adj.count(p))
                    for (int nb : pos_adj.at(p)) extended.insert(nb);
            }

            auto bfs_count = [&](const std::unordered_set<int>& nodes) -> int {
                std::unordered_set<int> vis; int comp = 0;
                for (int p : nodes) {
                    if (vis.count(p)) continue; comp++;
                    std::vector<int> q = {p}; vis.insert(p);
                    while (!q.empty()) {
                        int c = q.back(); q.pop_back();
                        if (pos_adj.count(c))
                            for (int nb : pos_adj.at(c))
                                if (nodes.count(nb) && !vis.count(nb))
                                    { vis.insert(nb); q.push_back(nb); }
                    }
                }
                return comp;
            };
            int before = bfs_count(extended);
            auto without = extended; without.erase(v.pos);
            if (!without.empty() && bfs_count(without) > before)
                cd.bridge_set.insert(v.pos);
        }

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

    int total = 0, unph = 0;
    for (auto& [ch, pm] : predictions_)
        for (auto& [pos, pr] : pm) {
            total++;
            if (pr.prob_error >= params_.break_threshold &&
                (!params_.respect_bridge || !pr.is_bridge)) unph++;
        }
    std::cerr << "\n[GNN] Predictions: " << total << ", unphase: " << unph << "\n";
    writeOutputVCF();
    if (!params_.sv_vcf.empty() && !params_.output_sv_vcf.empty())
        writeOutputVCF(params_.sv_vcf, params_.output_sv_vcf);
    if (!params_.mod_vcf.empty() && !params_.output_mod_vcf.empty())
        writeOutputVCF(params_.mod_vcf, params_.output_mod_vcf);
    return 0;
}

// ═══════════════════════════════════════════════════════════
void GNNModule::loadModel() {
    ort_options_.SetIntraOpNumThreads(1);
    ort_options_.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    ort_session_ = std::make_unique<Ort::Session>(
        ort_env_, params_.model_path.c_str(), ort_options_);
}

// ═══════════════════════════════════════════════════════════
void GNNModule::parseVCF() {
    htsFile* fp = hts_open(params_.vcf_path.c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
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
    for (auto& [ch, vl] : variants_)
        std::sort(vl.begin(), vl.end(),
                  [](const VariantInfo& a, const VariantInfo& b){ return a.pos < b.pos; });
    std::cerr << "  " << n << " phased variants\n";
}

// ═══════════════════════════════════════════════════════════
// DOT parsing: sscanf (was std::regex, 10-50x faster)
// ═══════════════════════════════════════════════════════════
void GNNModule::parseDotFiles() {
    // Collect chromosome list and pre-initialize dot_edges_ map
    std::vector<std::string> chroms;
    for (auto& [ch, vl] : variants_) {
        chroms.push_back(ch);
        dot_edges_[ch]; // pre-create entry to avoid concurrent map insertion
    }

    // Parse all chromosomes in parallel
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
                std::string path = params_.dot_prefix + "." + ch + ".dot";
                std::ifstream fin(path);
                if (!fin.is_open()) continue;

                std::vector<DotEdge> local_edges;
                std::string line;
                while (std::getline(fin, line)) {
                    int sp, sa, dp, da; float w;
                    if (std::sscanf(line.c_str(), "%d.%d -> %d.%d [label=%f]",
                                    &sp, &sa, &dp, &da, &w) == 5) {
                        local_edges.push_back({sp-1, sa, dp-1, da, w});
                    }
                }
                // Index by src_pos
                auto& em = dot_edges_[ch];
                for (auto& e : local_edges) em[e.src_pos].push_back(e);
                total += (int)local_edges.size();
            }
        });
    }
    for (auto& t : threads) t.join();
    std::cerr << "  " << total.load() << " DOT edges\n";
}

// ═══════════════════════════════════════════════════════════
void GNNModule::computeGenomicFeatures() {
    std::vector<std::string> chroms;
    for (auto& [ch, vl] : variants_) chroms.push_back(ch);

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
    if (N > MAX_NODES) return;

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

    // reverse edge weight
    std::unordered_map<int64_t, float> rev_w_map;
    auto rev_key = [](int dp, int da, int sp, int sa) -> int64_t {
        return ((int64_t)dp << 34) | ((int64_t)da << 32) | ((int64_t)sp << 4) | sa; };
    for (int pos : wpos_set) {
        if (!edge_map.count(pos)) continue;
        for (auto& e : edge_map.at(pos)) {
            if (!wpos_set.count(e.dst_pos)) continue;
            rev_w_map[rev_key(e.dst_pos, e.dst_allele, e.src_pos, e.src_allele)] = e.weight;
        }
    }

    // ── Node features [N × 27] ──
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
        float bridge_f = cd.bridge_set.count(v->pos) ? 1.0f : 0.0f;

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
            float rw = 0; auto rit = rev_w_map.find(rev_key(e.dst_pos,e.dst_allele,e.src_pos,e.src_allele));
            if (rit != rev_w_map.end()) rw = rit->second;
            int row = dn, col = sn;
            adj[row*N+col] = 1.0f;
            float* ep = &ef[(row*N+col)*EDGE_FEAT_DIM];
            ep[0] = 1.0f/(1.0f+std::exp(-e.weight));
            ep[1] = std::log10(1.0f+std::abs((float)(e.dst_pos-e.src_pos)));
            ep[2] = std::log1p(e.weight);
            ep[3] = same_hp; ep[4] = inter_ps;
            ep[5] = e.weight/(wsum+EPS); ep[6] = std::log1p(rw);
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
            for(int f=0;f<EDGE_FEAT_DIM;++f) se[f]+=ef[(i*N+j)*EDGE_FEAT_DIM+f]; ni++;
        }
        if (ni>0) for(int f=0;f<EDGE_FEAT_DIM;++f) ef[(i*N+i)*EDGE_FEAT_DIM+f]=se[f]/ni;
    }

    // ── ONNX inference ──
    std::array<int64_t,2> ns={N,NODE_FEAT_DIM}, as_={N,N};
    std::array<int64_t,3> es={N,N,EDGE_FEAT_DIM};
    auto tn = Ort::Value::CreateTensor<float>(memory_info_,nf.data(),nf.size(),ns.data(),2);
    auto ta = Ort::Value::CreateTensor<float>(memory_info_,adj.data(),adj.size(),as_.data(),2);
    auto te = Ort::Value::CreateTensor<float>(memory_info_,ef.data(),ef.size(),es.data(),3);
    std::vector<Ort::Value> in; in.reserve(3);
    in.push_back(std::move(tn)); in.push_back(std::move(ta)); in.push_back(std::move(te));
    const char* inn[]={"node_features","adjacency","edge_features"};
    const char* outn[]={"probabilities"};
    auto out = ort_session_->Run(Ort::RunOptions{nullptr},inn,in.data(),3,outn,1);
    float* op = out[0].GetTensorMutableData<float>();

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
        bool br = cd.bridge_set.count(wvars[i]->pos) > 0;
        std::lock_guard<std::mutex> lk(pred_mutex_);
        auto& pr = predictions_[chrom][wvars[i]->pos];
        pr.prob_error = (pr.prob_error*pr.count + p_err)/(pr.count+1);
        pr.is_bridge = pr.is_bridge || br;
        pr.count++;
    }
}

// ═══════════════════════════════════════════════════════════
void GNNModule::writeOutputVCF() {
    writeOutputVCF(params_.vcf_path, params_.output_vcf);
}

void GNNModule::writeOutputVCF(const std::string& in_path, const std::string& out_path) {
    htsFile* ifp = hts_open(in_path.c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(ifp);
    htsFile* ofp = hts_open(out_path.c_str(), "w");
    (void)bcf_hdr_write(ofp, hdr);
    bcf1_t* rec = bcf_init();
    int nu = 0;
    while (bcf_read(ifp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_ALL);
        std::string ch = bcf_hdr_id2name(hdr, rec->rid);
        bool unph = false;
        if (predictions_.count(ch) && predictions_[ch].count(rec->pos)) {
            auto& pr = predictions_[ch][rec->pos];
            if (pr.prob_error >= params_.break_threshold &&
                (!params_.respect_bridge || !pr.is_bridge)) unph = true;
        }
        if (unph) {
            int32_t* gt=nullptr; int ng=0;
            if (bcf_get_genotypes(hdr,rec,&gt,&ng)>=2) {
                gt[0]=bcf_gt_unphased(bcf_gt_allele(gt[0]));
                gt[1]=bcf_gt_unphased(bcf_gt_allele(gt[1]));
                bcf_update_genotypes(hdr,rec,gt,ng); free(gt);
            }
            bcf_update_format_int32(hdr,rec,"PS",nullptr,0);
            nu++;
        }
        (void)bcf_write(ofp, hdr, rec);
    }
    bcf_destroy(rec); bcf_hdr_destroy(hdr);
    hts_close(ifp); hts_close(ofp);
    std::cerr << "[GNN] " << out_path << ": unphased " << nu << " variants\n";
}

// ═══════════════════════════════════════════════════════════
void GNNModule::parseSecondaryVCF(const std::string& path, GnnVariantType vtype) {
    htsFile* fp = hts_open(path.c_str(), "r");
    if (!fp) { std::cerr << "  WARNING: cannot open " << path << "\n"; return; }
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
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
        vi.vtype = vtype;
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

// ═══════════════════════════════════════════════════════════
#ifdef GNN_STANDALONE
int main(int argc, char* argv[]) {
    GNNModule::Params p;
    for(int i=1;i<argc;++i){std::string a=argv[i];
        if((a=="--model"||a=="-m")&&i+1<argc) p.model_path=argv[++i];
        else if((a=="--vcf"||a=="-v")&&i+1<argc) p.vcf_path=argv[++i];
        else if(a=="--dot-prefix"&&i+1<argc) p.dot_prefix=argv[++i];
        else if((a=="--reference"||a=="-r")&&i+1<argc) p.reference_path=argv[++i];
        else if((a=="--output"||a=="-o")&&i+1<argc) p.output_vcf=argv[++i];
        else if((a=="-B"||a=="--break-threshold")&&i+1<argc) p.break_threshold=std::stof(argv[++i]);
        else if(a=="--pe-threshold"&&i+1<argc) p.pe_threshold=std::stof(argv[++i]);
        else if(a=="--window"&&i+1<argc) p.window=std::stoi(argv[++i]);
        else if((a=="-t"||a=="--threads")&&i+1<argc) p.threads=std::stoi(argv[++i]);
        else if(a=="--sv-vcf"&&i+1<argc) p.sv_vcf=argv[++i];
        else if(a=="--mod-vcf"&&i+1<argc) p.mod_vcf=argv[++i];
        else if(a=="--output-sv-vcf"&&i+1<argc) p.output_sv_vcf=argv[++i];
        else if(a=="--output-mod-vcf"&&i+1<argc) p.output_mod_vcf=argv[++i];
        else if(a=="--respect-bridge") p.respect_bridge=true;
        else if(a=="-h"||a=="--help"){std::cerr<<"longphase gnn --model M --vcf V --dot-prefix D -r R -o O [-B 0.30] [-t 4]\n"
            "  [--sv-vcf SV.vcf] [--mod-vcf MOD.vcf] [--output-sv-vcf O_SV] [--output-mod-vcf O_MOD]\n"
            "  [--respect-bridge]\n";return 0;}
    }
    if(p.model_path.empty()||p.vcf_path.empty()||p.output_vcf.empty()){std::cerr<<"Error: --model, --vcf, --output required\n";return 1;}
    return GNNModule(p).run();
}
#endif
