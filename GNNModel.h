// GNNModel.h -- inference for the phasing GPS model without ONNX Runtime.
//
// The weights are compiled in via GNNWeights.h, so no model file is needed at
// run time. The computation matches the exported ONNX graph. Everything is
// dense:
// the adjacency is [N, N] and the edge features [N, N, 7], so no sparse
// scatter primitives are required.
//
// Memory is O(N^2 * heads) for the attention scores, which is fine for the
// windows used here (tens to a couple of hundred nodes) but would need
// revisiting for whole-chromosome graphs.

#ifndef GNN_MODEL_H
#define GNN_MODEL_H

#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "GNNWeights.h"   // generated: dimensions + the parameter blob

namespace gnn {

// kNodeFeat, kEdgeFeat, kHidden, kHeads, kLayers, kClasses and kParamCount
// all come from the generated header.
constexpr int   kHeadDim  = kHidden / kHeads;
constexpr float kEps      = 1e-5f;      // LayerNorm and BatchNorm
constexpr float kLeaky    = 0.2f;       // GATv2 attention
constexpr float kMaskFill = -1e9f;      // masked-out attention logits
constexpr float kAdjThr   = 0.5f;       // adjacency values below this are absent

// One dense tensor, stored row-major.
struct Tensor {
    std::vector<float> data;
    int rows = 0, cols = 0;
    inline const float* row(int r) const { return data.data() + (size_t)r * cols; }
};

struct LayerWeights {
    Tensor lin_l_w, lin_r_w, edge_enc;      // [128,128] [128,128] [7,128]
    std::vector<float> lin_l_b, lin_r_b, att, mpnn_b;
    Tensor q_w, k_w, v_w, o_w;
    std::vector<float> q_b, k_b, v_b, o_b;
    std::vector<float> norm1_w, norm1_b, norm2_w, norm2_b;
    Tensor ffn0_w, ffn3_w;                  // [256,128] [128,256]
    std::vector<float> ffn0_b, ffn3_b;
};

class Model {
public:
    // Decodes the compiled-in weights. Cheap enough to do once at startup;
    // there is no file to find and nothing to go missing at deployment.
    Model() { load(); }

    // node_feat  [n * 31]     row-major
    // adjacency  [n * n]      adjacency(i,j) > 0.5 means an edge j -> i
    // edge_feat  [n * n * 7]  features of that same edge
    // returns    [n * 2]      class probabilities per node
    std::vector<float> forward(const std::vector<float>& node_feat,
                               const std::vector<float>& adjacency,
                               const std::vector<float>& edge_feat,
                               int n) const;

private:
    std::vector<float> in_w_, in_b_, in_mean_, in_var_;
    Tensor proj_w_;   std::vector<float> proj_b_;
    LayerWeights layer_[kLayers];
    Tensor cls0_w_, cls3_w_;
    std::vector<float> cls0_b_, cls3_b_;

    void load();
};

// ── elementwise helpers ─────────────────────────────────────────────────

inline float gelu(float x) {
    // Exact GELU, matching torch.nn.functional.gelu. The tanh approximation
    // differs by ~1e-3 and must not be substituted here.
    return 0.5f * x * (1.0f + std::erf(x * 0.70710678118654752f));
}

inline float leaky_relu(float x) { return x >= 0.0f ? x : kLeaky * x; }

// out[n, out_dim] = in[n, in_dim] * W^T + b, with W stored as [out_dim, in_dim]
inline void linear(const float* in, int n, int in_dim,
                   const Tensor& w, const std::vector<float>& b,
                   float* out) {
    const int out_dim = w.rows;
    for (int i = 0; i < n; ++i) {
        const float* xi = in + (size_t)i * in_dim;
        float* oi = out + (size_t)i * out_dim;
        for (int o = 0; o < out_dim; ++o) {
            const float* wo = w.row(o);
            float s = b[o];
            for (int k = 0; k < in_dim; ++k) s += xi[k] * wo[k];
            oi[o] = s;
        }
    }
}

inline void layer_norm(float* x, int n, int dim,
                       const std::vector<float>& w, const std::vector<float>& b) {
    for (int i = 0; i < n; ++i) {
        float* xi = x + (size_t)i * dim;
        float mean = 0.0f;
        for (int d = 0; d < dim; ++d) mean += xi[d];
        mean /= dim;
        float var = 0.0f;
        for (int d = 0; d < dim; ++d) { float t = xi[d] - mean; var += t * t; }
        var /= dim;                       // biased variance, as in PyTorch
        const float inv = 1.0f / std::sqrt(var + kEps);
        for (int d = 0; d < dim; ++d) xi[d] = (xi[d] - mean) * inv * w[d] + b[d];
    }
}

// Softmax along the second index of a [n, n, heads] score block.
inline void softmax_dim1(std::vector<float>& s, int n) {
    for (int i = 0; i < n; ++i) {
        for (int h = 0; h < kHeads; ++h) {
            float mx = -3.4e38f;
            for (int j = 0; j < n; ++j) {
                float v = s[((size_t)i * n + j) * kHeads + h];
                if (v > mx) mx = v;
            }
            float sum = 0.0f;
            for (int j = 0; j < n; ++j) {
                size_t idx = ((size_t)i * n + j) * kHeads + h;
                float e = std::exp(s[idx] - mx);
                s[idx] = e;
                sum += e;
            }
            const float inv = (sum > 0.0f) ? 1.0f / sum : 0.0f;
            for (int j = 0; j < n; ++j) s[((size_t)i * n + j) * kHeads + h] *= inv;
        }
    }
}

// ── implementation ──────────────────────────────────────────────────────

// Sequential reader over the decoded parameter blob. Running off the end
// means the header and this file disagree, which is a build error, so it
// throws rather than reading garbage.
namespace detail {
struct Cursor {
    const std::vector<float>& buf;
    size_t at = 0;
    explicit Cursor(const std::vector<float>& b) : buf(b) {}
    void take(Tensor& t, int rows, int cols) {
        const size_t n = (size_t)rows * cols;
        if (at + n > buf.size()) throw std::runtime_error("GNN weights exhausted");
        t.rows = rows; t.cols = cols;
        t.data.assign(buf.begin() + at, buf.begin() + at + n);
        at += n;
    }
    void take(std::vector<float>& v, int n) {
        if (at + (size_t)n > buf.size()) throw std::runtime_error("GNN weights exhausted");
        v.assign(buf.begin() + at, buf.begin() + at + n);
        at += n;
    }
};
}  // namespace detail

inline void Model::load() {
    const std::vector<float> blob = detail::decode_embedded_weights();
    if (blob.size() != kParamCount)
        throw std::runtime_error("embedded GNN weights are the wrong size");
    detail::Cursor c(blob);

    c.take(in_w_, kNodeFeat);
    c.take(in_b_, kNodeFeat);
    c.take(in_mean_, kNodeFeat);
    c.take(in_var_, kNodeFeat);
    c.take(proj_w_, kHidden, kNodeFeat);
    c.take(proj_b_, kHidden);

    for (int l = 0; l < kLayers; ++l) {
        LayerWeights& L = layer_[l];
        c.take(L.lin_l_w, kHidden, kHidden); c.take(L.lin_l_b, kHidden);
        c.take(L.lin_r_w, kHidden, kHidden); c.take(L.lin_r_b, kHidden);
        c.take(L.edge_enc, kEdgeFeat, kHidden);
        c.take(L.att, kHeads * kHeadDim);
        c.take(L.mpnn_b, kHidden);
        c.take(L.q_w, kHidden, kHidden); c.take(L.q_b, kHidden);
        c.take(L.k_w, kHidden, kHidden); c.take(L.k_b, kHidden);
        c.take(L.v_w, kHidden, kHidden); c.take(L.v_b, kHidden);
        c.take(L.o_w, kHidden, kHidden); c.take(L.o_b, kHidden);
        c.take(L.norm1_w, kHidden); c.take(L.norm1_b, kHidden);
        c.take(L.ffn0_w, kHidden * 2, kHidden); c.take(L.ffn0_b, kHidden * 2);
        c.take(L.ffn3_w, kHidden, kHidden * 2); c.take(L.ffn3_b, kHidden);
        c.take(L.norm2_w, kHidden); c.take(L.norm2_b, kHidden);
    }

    c.take(cls0_w_, kHidden, kHidden + 2); c.take(cls0_b_, kHidden);
    c.take(cls3_w_, kClasses, kHidden);    c.take(cls3_b_, kClasses);

    if (c.at != blob.size())
        throw std::runtime_error("embedded GNN weights not fully consumed");
}

inline std::vector<float> Model::forward(const std::vector<float>& node_feat,
                                         const std::vector<float>& adjacency,
                                         const std::vector<float>& edge_feat,
                                         int n) const {
    const size_t nn = (size_t)n * n;

    // Edge context: per node, the sum and the max of the raw edge weight over
    // incoming edges, self-loops excluded. Matches the ReduceSum / ReduceMax
    // pair at the head of the exported graph.
    std::vector<float> ctx((size_t)n * 2, 0.0f);
    for (int i = 0; i < n; ++i) {
        float sum = 0.0f, mx = kMaskFill;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            if (adjacency[(size_t)i * n + j] <= 0.0f) continue;
            const float w = edge_feat[((size_t)i * n + j) * kEdgeFeat];
            sum += w;
            if (w > mx) mx = w;
        }
        ctx[(size_t)i * 2 + 0] = sum;
        ctx[(size_t)i * 2 + 1] = mx > 0.0f ? mx : 0.0f;   // Clip(min=0)
    }

    // Input BatchNorm (inference: a fixed affine map) then projection + GELU.
    std::vector<float> normed((size_t)n * kNodeFeat);
    for (int i = 0; i < n; ++i)
        for (int d = 0; d < kNodeFeat; ++d)
            normed[(size_t)i * kNodeFeat + d] =
                (node_feat[(size_t)i * kNodeFeat + d] - in_mean_[d]) /
                std::sqrt(in_var_[d] + kEps) * in_w_[d] + in_b_[d];

    std::vector<float> x((size_t)n * kHidden);
    linear(normed.data(), n, kNodeFeat, proj_w_, proj_b_, x.data());
    for (float& v : x) v = gelu(v);

    // Scratch buffers reused across layers.
    // Row-major edge list: edges holds the flat index i*n+j of every present
    // edge, grouped by target row, and row_beg indexes into it.
    std::vector<size_t> edges;
    std::vector<int> row_beg(n + 1, 0);
    edges.reserve(nn / 4);
    for (int i = 0; i < n; ++i) {
        row_beg[i] = (int)edges.size();
        for (int j = 0; j < n; ++j)
            if (adjacency[(size_t)i * n + j] >= kAdjThr)
                edges.push_back((size_t)i * n + j);
    }
    row_beg[n] = (int)edges.size();

    std::vector<float> xl((size_t)n * kHidden), xr((size_t)n * kHidden);
    std::vector<float> edge_h(edges.size() * kHidden);
    std::vector<float> escore(edges.size() * kHeads);
    std::vector<float> scores(nn * kHeads);
    std::vector<float> local((size_t)n * kHidden), glob((size_t)n * kHidden);
    std::vector<float> q((size_t)n * kHidden), k((size_t)n * kHidden), v((size_t)n * kHidden);
    std::vector<float> agg((size_t)n * kHidden);
    std::vector<float> ffn_h((size_t)n * kHidden * 2);

    for (int l = 0; l < kLayers; ++l) {
        const LayerWeights& L = layer_[l];

        // ---- local branch: GATv2 over the dense adjacency ----
        linear(x.data(), n, kHidden, L.lin_l_w, L.lin_l_b, xl.data());
        linear(x.data(), n, kHidden, L.lin_r_w, L.lin_r_b, xr.data());

        // Only present edges are projected and scored. Absent entries would
        // receive kMaskFill and contribute a vanishing softmax weight, so
        // skipping them is exact as long as rows with no edge at all are
        // handled separately (see below).
        for (size_t t = 0; t < edges.size(); ++t) {
            const size_t e = edges[t];
            const float* ef = edge_feat.data() + e * kEdgeFeat;
            float* eh = edge_h.data() + (size_t)t * kHidden;
            for (int o = 0; o < kHidden; ++o) {
                float s = 0.0f;
                for (int c = 0; c < kEdgeFeat; ++c) s += ef[c] * L.edge_enc.row(c)[o];
                eh[o] = s;
            }
        }

        // score(i,j,h) = att . leaky_relu(x_r[i] + x_l[j] + edge[i,j])
        // Row i is the target, column j the source.
        std::fill(agg.begin(), agg.end(), 0.0f);
        for (int i = 0; i < n; ++i) {
            float* ai = agg.data() + (size_t)i * kHidden;
            const int beg = row_beg[i], end = row_beg[i + 1];

            if (beg == end) {
                // No incoming edge: every logit is kMaskFill, so the softmax
                // is uniform over all n columns. Reproduce that exactly.
                const float w = 1.0f / (float)n;
                for (int j = 0; j < n; ++j) {
                    const float* xlj = xl.data() + (size_t)j * kHidden;
                    for (int d = 0; d < kHidden; ++d) ai[d] += w * xlj[d];
                }
                for (int d = 0; d < kHidden; ++d) ai[d] += L.mpnn_b[d];
                continue;
            }

            const float* xri = xr.data() + (size_t)i * kHidden;
            for (int t = beg; t < end; ++t) {
                const int j = (int)(edges[t] % (size_t)n);
                const float* xlj = xl.data() + (size_t)j * kHidden;
                const float* eh  = edge_h.data() + (size_t)t * kHidden;
                for (int h = 0; h < kHeads; ++h) {
                    const int off = h * kHeadDim;
                    float s = 0.0f;
                    for (int d = 0; d < kHeadDim; ++d)
                        s += leaky_relu(xri[off + d] + xlj[off + d] + eh[off + d]) *
                             L.att[off + d];
                    escore[(size_t)t * kHeads + h] = s;
                }
            }
            // softmax per head over this row's edges
            for (int h = 0; h < kHeads; ++h) {
                float mx = -3.4e38f;
                for (int t = beg; t < end; ++t)
                    mx = std::max(mx, escore[(size_t)t * kHeads + h]);
                float sum = 0.0f;
                for (int t = beg; t < end; ++t) {
                    float e2 = std::exp(escore[(size_t)t * kHeads + h] - mx);
                    escore[(size_t)t * kHeads + h] = e2;
                    sum += e2;
                }
                const float inv = (sum > 0.0f) ? 1.0f / sum : 0.0f;
                const int off = h * kHeadDim;
                for (int t = beg; t < end; ++t) {
                    const float w = escore[(size_t)t * kHeads + h] * inv;
                    const float* xlj = xl.data() + (edges[t] % (size_t)n) * kHidden;
                    for (int d = 0; d < kHeadDim; ++d) ai[off + d] += w * xlj[off + d];
                }
            }
            for (int d = 0; d < kHidden; ++d) ai[d] += L.mpnn_b[d];
        }
        local = agg;

        // ---- global branch: multi-head self-attention over the window ----
        linear(x.data(), n, kHidden, L.q_w, L.q_b, q.data());
        linear(x.data(), n, kHidden, L.k_w, L.k_b, k.data());
        linear(x.data(), n, kHidden, L.v_w, L.v_b, v.data());

        const float scale = 1.0f / std::sqrt((float)kHeadDim);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int h = 0; h < kHeads; ++h) {
                    const int off = h * kHeadDim;
                    float s = 0.0f;
                    for (int d = 0; d < kHeadDim; ++d)
                        s += q[(size_t)i * kHidden + off + d] *
                             k[(size_t)j * kHidden + off + d];
                    scores[((size_t)i * n + j) * kHeads + h] = s * scale;
                }
        softmax_dim1(scores, n);

        std::fill(agg.begin(), agg.end(), 0.0f);
        for (int i = 0; i < n; ++i) {
            float* ai = agg.data() + (size_t)i * kHidden;
            for (int j = 0; j < n; ++j) {
                const float* vj = v.data() + (size_t)j * kHidden;
                const float* a = &scores[((size_t)i * n + j) * kHeads];
                for (int h = 0; h < kHeads; ++h) {
                    const float w = a[h];
                    const int off = h * kHeadDim;
                    for (int d = 0; d < kHeadDim; ++d) ai[off + d] += w * vj[off + d];
                }
            }
        }
        linear(agg.data(), n, kHidden, L.o_w, L.o_b, glob.data());

        // ---- parallel combine, then FFN ----
        for (size_t t = 0; t < x.size(); ++t) x[t] += local[t] + glob[t];
        layer_norm(x.data(), n, kHidden, L.norm1_w, L.norm1_b);

        linear(x.data(), n, kHidden, L.ffn0_w, L.ffn0_b, ffn_h.data());
        for (float& t : ffn_h) t = gelu(t);
        linear(ffn_h.data(), n, kHidden * 2, L.ffn3_w, L.ffn3_b, agg.data());
        for (size_t t = 0; t < x.size(); ++t) x[t] += agg[t];
        layer_norm(x.data(), n, kHidden, L.norm2_w, L.norm2_b);
    }

    // ---- classifier over [hidden | edge context] ----
    std::vector<float> cat((size_t)n * (kHidden + 2));
    for (int i = 0; i < n; ++i) {
        std::memcpy(&cat[(size_t)i * (kHidden + 2)], &x[(size_t)i * kHidden],
                    kHidden * sizeof(float));
        cat[(size_t)i * (kHidden + 2) + kHidden]     = ctx[(size_t)i * 2];
        cat[(size_t)i * (kHidden + 2) + kHidden + 1] = ctx[(size_t)i * 2 + 1];
    }
    std::vector<float> h1((size_t)n * kHidden);
    linear(cat.data(), n, kHidden + 2, cls0_w_, cls0_b_, h1.data());
    for (float& t : h1) t = gelu(t);

    std::vector<float> logits((size_t)n * kClasses);
    linear(h1.data(), n, kHidden, cls3_w_, cls3_b_, logits.data());

    for (int i = 0; i < n; ++i) {
        float* li = logits.data() + (size_t)i * kClasses;
        float mx = li[0];
        for (int c = 1; c < kClasses; ++c) mx = std::max(mx, li[c]);
        float sum = 0.0f;
        for (int c = 0; c < kClasses; ++c) { li[c] = std::exp(li[c] - mx); sum += li[c]; }
        for (int c = 0; c < kClasses; ++c) li[c] /= sum;
    }
    return logits;
}

}  // namespace gnn

#endif  // GNN_MODEL_H
