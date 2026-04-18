// =====================================================================
//  CompareProcess.cpp
//  -------------------------------------------------------------------
//  LongPhase "compare" sub-command – processing core.
//
//  Adapted by Claude (Anthropic) from WhatsHap's compare module:
//      https://github.com/whatshap/whatshap/blob/main/whatshap/cli/compare.py
//
//  Implementation notes
//  --------------------
//  * VCF parsing is done via htslib (already a dependency of longphase).
//  * For every chromosome we build the intersection of heterozygous
//    phased variants of the truth (VCF 0) and the query (VCF 1).
//  * Variants in the intersection inherit their phase-set ID from each
//    VCF; blocks are identified by pairing the two PS values, i.e. the
//    intersection is taken over (ps_truth, ps_query) tuples
//    (same as WhatsHap).
//  * For every intersection block we compute in parallel:
//        - switch errors               (hamming of switch-encoded strings)
//        - SNV-only switch errors      (same but after filtering SNVs)
//        - minimum block-wise Hamming  (min over the two possible hap
//                                       assignments, diploid)
//        - switch / flip decomposition (WhatsHap compute_switch_flips)
//  * Block span statistics are computed on the QUERY VCF, chromosome by
//    chromosome.  N50 is over block spans in bp.
//  * At the bottom of the TSV we emit the two summary lines requested
//    by the user.
// =====================================================================

#include "CompareProcess.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <htslib/hts.h>
#include <htslib/vcf.h>

// ====================================================================
//  Small helpers
// ====================================================================

// a|b  with ploidy 2 – convert to a 2-char phase string, or "" if un-phased
static inline std::string percentStr(double num, double den)
{
    if (den <= 0) return "--";
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.5f", num * 100.0 / den);
    return buf;
}

// Strip directory prefix from a file path
static inline std::string baseName(const std::string& path)
{
    auto pos = path.find_last_of("/\\");
    return (pos == std::string::npos) ? path : path.substr(pos + 1);
}

// "switch encoding": given a bit string s, produce a bit string of
// length len(s)-1 where position i is 1 iff s[i] != s[i+1].
static std::string switchEncoding(const std::string& s)
{
    if (s.size() < 2) return "";
    std::string out(s.size() - 1, '0');
    for (size_t i = 1; i < s.size(); ++i)
        out[i - 1] = (s[i - 1] == s[i]) ? '0' : '1';
    return out;
}

// Hamming distance between two equal-length strings.
static inline int hammingStr(const std::string& a, const std::string& b)
{
    int d = 0;
    const size_t n = std::min(a.size(), b.size());
    for (size_t i = 0; i < n; ++i) if (a[i] != b[i]) ++d;
    return d;
}

// WhatsHap's compute_switch_flips (diploid, hap 0 only).
// Returns {switches, flips}.
static std::pair<int,int> computeSwitchFlips(const std::string& p0,
                                             const std::string& p1)
{
    std::string s0 = switchEncoding(p0);
    std::string s1 = switchEncoding(p1);
    int switches = 0, flips = 0;
    int switches_in_a_row = 0;
    const size_t n = s0.size();
    for (size_t i = 0; i < n; ++i) {
        if (s0[i] != s1[i]) ++switches_in_a_row;
        if (i + 1 == n || s0[i] == s1[i]) {
            flips    += switches_in_a_row / 2;
            switches += switches_in_a_row % 2;
            switches_in_a_row = 0;
        }
    }
    return {switches, flips};
}

// ====================================================================
//  VCF loading – using htslib
// ====================================================================

bool CompareProcess::loadVcf(const std::string& path,
                             std::map<std::string, std::vector<VariantRecord>>& out,
                             std::string& sampleUsed,
                             long long& outTotalPhasedSnv)
{
    outTotalPhasedSnv = 0;

    htsFile* fp = bcf_open(path.c_str(), "r");
    if (!fp) {
        std::cerr << "[compare] ERROR: cannot open " << path << "\n";
        return false;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        std::cerr << "[compare] ERROR: cannot read header of " << path << "\n";
        bcf_close(fp);
        return false;
    }

    // ---- pick sample index ----
    int nSamples = bcf_hdr_nsamples(hdr);
    if (nSamples < 1) {
        std::cerr << "[compare] ERROR: no samples in " << path << "\n";
        bcf_hdr_destroy(hdr); bcf_close(fp);
        return false;
    }
    int sampleIdx = 0;
    if (!params.sample.empty()) {
        int idx = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, params.sample.c_str());
        if (idx < 0) {
            if (params.ignoreSampleName) {
                sampleIdx = 0;
            } else {
                std::cerr << "[compare] ERROR: sample '" << params.sample
                          << "' not found in " << path << "\n";
                bcf_hdr_destroy(hdr); bcf_close(fp);
                return false;
            }
        } else sampleIdx = idx;
    }
    sampleUsed = hdr->samples[sampleIdx];

    // ---- iterate records ----
    bcf1_t* rec = bcf_init();
    int32_t* gt_arr = nullptr;  int gt_n = 0;
    int32_t* ps_arr = nullptr;  int ps_n = 0;

    // Counter: how many phased variants lacked a PS tag and therefore
    // received a synthetic block id (= 1, i.e. whole-chromosome block).
    // Common for GIAB / T2T benchmark VCFs.
    size_t synthPsCount = 0;
    bool   psFieldSeen  = false;   // whether the header even declares PS

    while (bcf_read(fp, hdr, rec) >= 0) {
        bcf_unpack(rec, BCF_UN_ALL);

        // Must have REF + at least one ALT
        if (rec->n_allele < 2) continue;

        // We only look at the first ALT (diploid bi-allelic focus, as
        // WhatsHap does in the default configuration).
        const char* ref = rec->d.allele[0];
        const char* alt = rec->d.allele[1];
        VariantRecord v;
        v.pos   = rec->pos;            // 0-based
        v.ref   = ref;
        v.alt   = alt;
        v.isSnv = (std::strlen(ref) == 1 && std::strlen(alt) == 1);

        if (params.onlySnvs && !v.isSnv) continue;

        // ---- GT ----
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &gt_n);
        if (ngt < 0) continue;
        int ploidy = ngt / nSamples;
        if (ploidy != 2) continue;          // strictly diploid
        int32_t* g = gt_arr + sampleIdx * ploidy;
        if (bcf_gt_is_missing(g[0]) || bcf_gt_is_missing(g[1])) continue;
        int a0 = bcf_gt_allele(g[0]);
        int a1 = bcf_gt_allele(g[1]);
        // Must be heterozygous between REF(0) and ALT(1); skip everything else
        if (!((a0 == 0 && a1 == 1) || (a0 == 1 && a1 == 0))) continue;
        // Phased?
        bool phased = (bcf_gt_is_phased(g[1]) != 0);
        v.h0 = phased ? (int8_t)a0 : (int8_t)-1;
        v.h1 = phased ? (int8_t)a1 : (int8_t)-1;

        // ---- PS ----
        v.ps = 0;
        if (phased) {
            int nps = bcf_get_format_int32(hdr, rec, "PS", &ps_arr, &ps_n);
            if (nps > 0) {
                psFieldSeen = true;
                int32_t pv = ps_arr[sampleIdx];
                if (pv != bcf_int32_missing && pv != 0) v.ps = pv;
            }
            // GT is phased but PS is absent (or 0 / missing) ⇒ treat the
            // whole chromosome as a single phase block.  This is the
            // standard convention for GIAB / T2T benchmark VCFs.
            if (v.ps == 0) {
                v.ps = 1;
                ++synthPsCount;
            }
        }

        // Tally every phased het SNV in this VCF (regardless of whether
        // it will be intersected later with another VCF).
        if (phased && v.isSnv) ++outTotalPhasedSnv;

        const char* chrom = bcf_hdr_id2name(hdr, rec->rid);
        out[chrom].push_back(v);
    }

    if (synthPsCount > 0) {
        std::cerr << "[compare] Notice: " << synthPsCount
                  << " phased variants in " << path
                  << (psFieldSeen
                      ? " had no PS value"
                      : " had no PS tag")
                  << "; treating each chromosome as a single phase block "
                     "(standard behaviour for benchmark VCFs).\n";
    }

    if (gt_arr) free(gt_arr);
    if (ps_arr) free(ps_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    return true;
}

// ====================================================================
//  N50 over block spans
// ====================================================================
long long CompareProcess::computeN50(std::vector<long long> spans)
{
    if (spans.empty()) return 0;
    std::sort(spans.begin(), spans.end(), std::greater<long long>());
    long long total = 0;
    for (auto s : spans) total += s;
    long long acc = 0;
    for (auto s : spans) {
        acc += s;
        if (acc * 2 >= total) return s;
    }
    return spans.back();
}

// ====================================================================
//  Per-chromosome comparison (multi-threaded over blocks)
// ====================================================================
ChromResult CompareProcess::compareChromosome(
        const std::string& chrom,
        const std::vector<VariantRecord>& truth,
        const std::vector<VariantRecord>& query)
{
    ChromResult R;  R.chromosome = chrom;

    // -------- 1. build look-ups keyed by (pos, ref, alt) ----------- //
    auto makeKey = [](const VariantRecord& v) {
        return VariantKey{v.pos, v.ref, v.alt};
    };
    std::map<VariantKey, const VariantRecord*> truthMap, queryMap;
    for (const auto& v : truth) truthMap[makeKey(v)] = &v;
    for (const auto& v : query) queryMap[makeKey(v)] = &v;

    // -------- 2. simple per-kind counts on the TRUTH set ----------- //
    for (const auto& v : truth) {
        if (v.isSnv) ++R.totalSnvTruth;
        else         ++R.totalIndelTruth;
    }

    // -------- 3. variants common to both VCFs ---------------------- //
    struct CommonVar { const VariantRecord* t; const VariantRecord* q; };
    std::vector<CommonVar> common;
    common.reserve(std::min(truth.size(), query.size()));
    for (const auto& v : truth) {
        auto it = queryMap.find(makeKey(v));
        if (it == queryMap.end()) continue;
        common.push_back({&v, it->second});
        ++R.commonHet;
        bool phasedInBoth = (v.h0 >= 0) && (it->second->h0 >= 0);
        if (phasedInBoth) {
            if (v.isSnv) ++R.phasedSnv;
            else         ++R.phasedIndel;
        }
    }

    // Sort common by position so block strings follow genome order.
    std::sort(common.begin(), common.end(),
              [](const CommonVar& a, const CommonVar& b) {
                  return a.t->pos < b.t->pos;
              });

    // -------- 4. block statistics on the QUERY ------------------- //
    {
        std::unordered_map<int32_t, std::vector<int32_t>> psPos;
        for (const auto& v : query) {
            if (v.ps == 0 || v.h0 < 0) continue;
            psPos[v.ps].push_back(v.pos);
        }
        std::vector<long long> spans;
        spans.reserve(psPos.size());
        for (auto& kv : psPos) {
            auto& vec = kv.second;
            if (vec.size() < 2) continue;
            auto mm = std::minmax_element(vec.begin(), vec.end());
            long long span = (long long)(*mm.second - *mm.first);
            spans.push_back(span);
            R.blockSum += span;
        }
        R.numBlocks  = (int)spans.size();
        R.blockN50   = computeN50(spans);
        R.blockSpans = std::move(spans);   // keep for global N50
    }

    // -------- 5. build intersection blocks on phased-in-both vars -- //
    // key: (ps_truth, ps_query) -> list of indices into 'common'
    struct PairPS {
        int32_t t; int32_t q;
        bool operator==(const PairPS& o) const { return t == o.t && q == o.q; }
    };
    struct PairPSHash {
        size_t operator()(const PairPS& p) const noexcept {
            return std::hash<uint64_t>()((uint64_t)(uint32_t)p.t << 32
                                         | (uint32_t)p.q);
        }
    };
    std::unordered_map<PairPS, std::vector<int>, PairPSHash> blocks;
    for (size_t i = 0; i < common.size(); ++i) {
        const auto& cv = common[i];
        if (cv.t->h0 < 0 || cv.q->h0 < 0) continue;
        if (cv.t->ps == 0 || cv.q->ps == 0) continue;
        blocks[{cv.t->ps, cv.q->ps}].push_back((int)i);
    }

    // Materialise the non-singleton blocks for the parallel workers.
    std::vector<std::vector<int>> workBlocks;
    workBlocks.reserve(blocks.size());
    for (auto& kv : blocks) {
        if (kv.second.size() >= 2) workBlocks.emplace_back(std::move(kv.second));
    }
    R.intersectionBlocks = (int)workBlocks.size();
    for (auto& b : workBlocks) R.coveredVariants += (int)b.size();

    // -------- 6. parallel switch / Hamming over blocks ------------- //
    struct BlockStat {
        int switches       = 0;
        int snvSwitches    = 0;
        int snvAssessedPairs = 0;
        int hamming        = 0;
        int compared       = 0;
        int assessedPairs  = 0;
        int sfSwitches     = 0;
        int sfFlips        = 0;
        std::vector<SwBedRecord> bedRecords;    // only filled if requested
    };
    std::vector<BlockStat> perBlock(workBlocks.size());

    const bool wantBed = !params.swBedFile.empty();

    int nThreads = std::max(1, params.numThreads);
    if ((size_t)nThreads > workBlocks.size())
        nThreads = (int)std::max<size_t>(1, workBlocks.size());

    std::atomic<size_t> nextIdx{0};
    auto worker = [&]() {
        while (true) {
            size_t i = nextIdx.fetch_add(1);
            if (i >= workBlocks.size()) break;
            const auto& block = workBlocks[i];
            BlockStat& bs = perBlock[i];

            // Build h0 / h1 strings for truth and query.
            std::string t0, t1, q0, q1;
            t0.reserve(block.size()); t1.reserve(block.size());
            q0.reserve(block.size()); q1.reserve(block.size());
            std::vector<char> isSnv; isSnv.reserve(block.size());
            for (int idx : block) {
                const auto& cv = common[idx];
                t0.push_back('0' + cv.t->h0);
                t1.push_back('0' + cv.t->h1);
                q0.push_back('0' + cv.q->h0);
                q1.push_back('0' + cv.q->h1);
                isSnv.push_back(cv.t->isSnv ? 1 : 0);
            }

            // -- Hamming: take min over the two haplotype assignments --
            int dSame = hammingStr(t0, q0) + hammingStr(t1, q1);
            int dFlip = hammingStr(t0, q1) + hammingStr(t1, q0);
            bs.hamming       = std::min(dSame, dFlip) / 2;
            bs.compared      = (int)block.size();
            bs.assessedPairs = (int)block.size() - 1;

            // -- switch error on h0 (identical for h1 in diploid) --
            std::string sT = switchEncoding(t0);
            std::string sQ = switchEncoding(q0);
            bs.switches = hammingStr(sT, sQ);

            // -- switch/flip decomposition --
            auto sf = computeSwitchFlips(t0, q0);
            bs.sfSwitches = sf.first;
            bs.sfFlips    = sf.second;

            // -- SNV-only switch: rebuild strings restricted to SNV --
            std::string t0s, q0s;
            t0s.reserve(block.size()); q0s.reserve(block.size());
            for (size_t k = 0; k < block.size(); ++k) {
                if (isSnv[k]) { t0s.push_back(t0[k]); q0s.push_back(q0[k]); }
            }
            if (t0s.size() >= 2) {
                bs.snvSwitches       = hammingStr(switchEncoding(t0s),
                                                  switchEncoding(q0s));
                bs.snvAssessedPairs  = (int)t0s.size() - 1;
            }

            // -- BED records: emit one per position where the two
            //    switch encodings disagree (adapted from WhatsHap
            //    BedCreator.records) ------------------------------ //
            if (wantBed) {
                for (size_t k = 0; k < sT.size(); ++k) {
                    if (sT[k] == sQ[k]) continue;
                    const auto& cvA = common[block[k]];       // position i
                    const auto& cvB = common[block[k + 1]];   // position i+1
                    SwBedRecord r;
                    r.chromosome = chrom;
                    r.start      = cvA.t->pos;                 // 0-based
                    r.end        = cvB.t->pos;                 // exclusive
                    r.isSnv      = cvA.t->isSnv && cvB.t->isSnv;
                    r.psTruth    = cvA.t->ps;
                    r.psQuery    = cvA.q->ps;
                    bs.bedRecords.push_back(std::move(r));
                }
            }
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(nThreads);
    for (int i = 0; i < nThreads; ++i) threads.emplace_back(worker);
    for (auto& t : threads) t.join();

    // -------- 7. accumulate per-block stats into ChromResult ------- //
    for (auto& bs : perBlock) {
        R.allSwitches        += bs.switches;
        R.assessedPairs      += bs.assessedPairs;
        R.snvSwitches        += bs.snvSwitches;
        R.snvAssessedPairs   += bs.snvAssessedPairs;
        R.blockwiseHamming   += bs.hamming;
        R.totalCompared      += bs.compared;
        R.switchFlipsSwitches += bs.sfSwitches;
        R.switchFlipsFlips    += bs.sfFlips;
        if (!bs.bedRecords.empty()) {
            R.swBedRecords.insert(R.swBedRecords.end(),
                std::make_move_iterator(bs.bedRecords.begin()),
                std::make_move_iterator(bs.bedRecords.end()));
        }
    }

    return R;
}

// ====================================================================
//  Report writing
// ====================================================================
void CompareProcess::writeReport(const std::vector<ChromResult>& chromResults,
                                 const std::string& outPath,
                                 const std::string& queryBasename,
                                 long long queryTotalPhasedSnv,
                                 long long queryTotalBlocks,
                                 long long queryTotalBlockN50,
                                 long long queryTotalBlockSum)
{
    std::ofstream ofs(outPath);
    if (!ofs) {
        std::cerr << "[compare] ERROR: cannot write " << outPath << "\n";
        return;
    }

    // -------------------------------------------------------------
    // 1.  Accumulate intersection-based totals (sum over per-chrom)
    // -------------------------------------------------------------
    long long tot_phasedSnv = 0, tot_totalSnv = 0;
    long long tot_phasedIndel = 0, tot_totalIndel = 0;
    long long tot_snvSw = 0, tot_snvSwPairs = 0;
    long long tot_hamming = 0, tot_compared = 0;
    long long tot_blocks = 0, tot_blockSum = 0;
    std::vector<long long> allBlockSpans;

    for (const auto& R : chromResults) {
        tot_phasedSnv   += R.phasedSnv;
        tot_totalSnv    += R.totalSnvTruth;
        tot_phasedIndel += R.phasedIndel;
        tot_totalIndel  += R.totalIndelTruth;
        tot_snvSw       += R.snvSwitches;
        tot_snvSwPairs  += R.snvAssessedPairs;
        tot_hamming     += R.blockwiseHamming;
        tot_compared    += R.totalCompared;
        tot_blocks      += R.numBlocks;
        tot_blockSum    += R.blockSum;
        allBlockSpans.insert(allBlockSpans.end(),
                             R.blockSpans.begin(), R.blockSpans.end());
    }
    long long globalN50 = computeN50(allBlockSpans);

    // -------------------------------------------------------------
    // 2.  Top-of-file metadata: raw totals from the QUERY VCF itself
    //     (independent of the intersection with the truth VCF).
    // -------------------------------------------------------------
    ofs << "##Total phased SNV = " << queryTotalPhasedSnv << "\n";
    ofs << "##Total Block = "      << queryTotalBlocks    << "\n";
    ofs << "##Total Block N50 = "  << queryTotalBlockN50  << "\n";
    ofs << "##Total Block Sum = "  << queryTotalBlockSum  << "\n";

    // -------------------------------------------------------------
    // 3.  ### Summary header + one data row (across all common chroms)
    // -------------------------------------------------------------
    ofs << "###Sample"
        << "\tPhased_SNV\tPhased_SNV(%)"
        << "\tPhased_INDEL\tPhased_INDEL(%)"
        << "\tSNV_SW\tSNV_SW(%)"
        << "\tHamming_Dis.(%)"
        << "\tNo._of_Block\tBlock_N50(bp)\tBlock_Sum"
        << "\n";

    ofs << "###" << queryBasename
        << "\t" << tot_phasedSnv
        << "\t" << percentStr(tot_phasedSnv, tot_totalSnv)
        << "\t" << tot_phasedIndel
        << "\t" << percentStr(tot_phasedIndel, tot_totalIndel)
        << "\t" << tot_snvSw
        << "\t" << percentStr(tot_snvSw, tot_snvSwPairs)
        << "\t" << percentStr(tot_hamming, tot_compared)
        << "\t" << tot_blocks
        << "\t" << globalN50
        << "\t" << tot_blockSum
        << "\n";

    // -------------------------------------------------------------
    // 4.  Per-chromosome header + one row per chromosome
    // -------------------------------------------------------------
    ofs << "#sample\tchromosome\tdataset1\tdataset2"
        << "\tcommon_het"
        << "\tphased_SNV\tphased_SNV(%)"
        << "\tphased_INDEL\tphased_INDEL(%)"
        << "\tintersection_blocks\tcovered_variants"
        << "\tall_assessed_pairs\tall_switches\tall_switch_rate(%)"
        << "\tswitchflips(sw/fl)\tswitchflip_rate(%)"
        << "\tblockwise_hamming\tblockwise_hamming_rate(%)"
        << "\tSNV_SW\tSNV_SW(%)"
        << "\tNo._of_Block\tBlock_N50(bp)\tBlock_Sum"
        << "\n";

    const std::string& d1 = params.datasetNames[0];
    const std::string& d2 = params.datasetNames[1];

    for (const auto& R : chromResults) {
        ofs << queryBasename << '\t' << R.chromosome << '\t' << d1 << '\t' << d2
            << '\t' << R.commonHet
            << '\t' << R.phasedSnv
            << '\t' << percentStr(R.phasedSnv,  R.totalSnvTruth)
            << '\t' << R.phasedIndel
            << '\t' << percentStr(R.phasedIndel, R.totalIndelTruth)
            << '\t' << R.intersectionBlocks
            << '\t' << R.coveredVariants
            << '\t' << R.assessedPairs
            << '\t' << R.allSwitches
            << '\t' << percentStr(R.allSwitches, R.assessedPairs)
            << '\t' << R.switchFlipsSwitches << '/' << R.switchFlipsFlips
            << '\t' << percentStr(R.switchFlipsSwitches + R.switchFlipsFlips,
                                  R.assessedPairs)
            << '\t' << R.blockwiseHamming
            << '\t' << percentStr(R.blockwiseHamming, R.totalCompared)
            << '\t' << R.snvSwitches
            << '\t' << percentStr(R.snvSwitches, R.snvAssessedPairs)
            << '\t' << R.numBlocks
            << '\t' << R.blockN50
            << '\t' << R.blockSum
            << '\n';
    }

    ofs.close();
    std::cerr << "[compare] report written to " << outPath << "\n";
}

// ====================================================================
//  run – top-level driver
// ====================================================================
int CompareProcess::run()
{
    // Currently support pairwise comparison (truth vs query).
    if (params.vcfFiles.size() != 2) {
        std::cerr << "[compare] Notice: " << params.vcfFiles.size()
                  << " VCFs were given; only the first two will be compared "
                     "(multiway comparison is not implemented).\n";
    }
    const std::string& truthPath = params.vcfFiles[0];
    const std::string& queryPath = params.vcfFiles[1];

    std::map<std::string, std::vector<VariantRecord>> truth, query;
    std::string sTruth, sQuery;
    long long truthTotalPhasedSnv = 0, queryTotalPhasedSnv = 0;

    std::cerr << "[compare] reading truth VCF : " << truthPath << "\n";
    if (!loadVcf(truthPath, truth, sTruth, truthTotalPhasedSnv)) return 1;

    std::cerr << "[compare] reading query VCF : " << queryPath << "\n";
    if (!loadVcf(queryPath, query, sQuery, queryTotalPhasedSnv)) return 1;

    if (!params.ignoreSampleName && sTruth != sQuery) {
        std::cerr << "[compare] Warning: sample names differ ("
                  << sTruth << " vs " << sQuery << ").  Use --ignore-sample-name "
                     "to suppress this warning.\n";
    }

    // ---- compute block totals over the ENTIRE query VCF --------- //
    //      (independent of which chromosomes intersect the truth)
    long long queryTotalBlocks = 0, queryTotalBlockSum = 0, queryTotalBlockN50 = 0;
    {
        std::vector<long long> allSpans;
        for (const auto& kv : query) {
            std::unordered_map<int32_t, std::vector<int32_t>> psPos;
            for (const auto& v : kv.second) {
                if (v.h0 < 0 || v.ps == 0) continue;
                psPos[v.ps].push_back(v.pos);
            }
            for (auto& p : psPos) {
                if (p.second.size() < 2) continue;
                auto mm = std::minmax_element(p.second.begin(), p.second.end());
                long long span = (long long)(*mm.second - *mm.first);
                allSpans.push_back(span);
                queryTotalBlockSum += span;
                ++queryTotalBlocks;
            }
        }
        queryTotalBlockN50 = computeN50(allSpans);
    }

    // intersection of chromosomes
    std::vector<std::string> chroms;
    for (auto& kv : truth) if (query.count(kv.first)) chroms.push_back(kv.first);
    std::sort(chroms.begin(), chroms.end());

    std::cerr << "[compare] " << chroms.size()
              << " common chromosomes; using " << params.numThreads
              << " thread(s).\n";

    std::vector<ChromResult> results;
    results.reserve(chroms.size());
    for (const auto& c : chroms) {
        std::cerr << "[compare] -- chromosome " << c << " --\n";
        results.push_back(compareChromosome(c, truth[c], query[c]));
    }

    writeReport(results,
                params.outPrefix + ".tsv",
                baseName(queryPath),
                queryTotalPhasedSnv,
                queryTotalBlocks,
                queryTotalBlockN50,
                queryTotalBlockSum);

    // ---- switch-error BED output (merge + sort + write) ---------- //
    if (!params.swBedFile.empty()) {
        std::vector<SwBedRecord> all;
        size_t total = 0;
        for (const auto& R : results) total += R.swBedRecords.size();
        all.reserve(total);
        for (auto& R : results) {
            all.insert(all.end(),
                       std::make_move_iterator(R.swBedRecords.begin()),
                       std::make_move_iterator(R.swBedRecords.end()));
        }
        std::sort(all.begin(), all.end());

        std::ofstream bed(params.swBedFile);
        if (!bed) {
            std::cerr << "[compare] ERROR: cannot write "
                      << params.swBedFile << "\n";
            return 1;
        }
        bed << "#chrom\tstart\tend\ttype\tps_truth\tps_query\tdataset_pair\n";
        const std::string pair = params.datasetNames[0] + "<->" + params.datasetNames[1];
        for (const auto& r : all) {
            bed << r.chromosome << '\t' << r.start << '\t' << r.end
                << '\t' << (r.isSnv ? "SNV" : "INDEL")
                << '\t' << r.psTruth << '\t' << r.psQuery
                << '\t' << pair << '\n';
        }
        std::cerr << "[compare] BED with " << all.size()
                  << " switch errors written to "
                  << params.swBedFile << "\n";
    }

    return 0;
}
