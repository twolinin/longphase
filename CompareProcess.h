// =====================================================================
//  CompareProcess.h
//  -------------------------------------------------------------------
//  LongPhase "compare" sub-command – processing core.
//
//  Adapted by Claude (Anthropic) from WhatsHap's compare module:
//      https://github.com/whatshap/whatshap/blob/main/whatshap/cli/compare.py
//
//  This header declares:
//    - CompareParameters  : holds all CLI parameters
//    - VariantRecord      : one het variant with phasing & PS info
//    - ChromResult        : per-chromosome comparison result
//    - CompareProcess     : driver class
// =====================================================================
#ifndef COMPAREPROCESS_H
#define COMPAREPROCESS_H

#include <cstdint>
#include <map>
#include <string>
#include <vector>

// ------------------------------------------------------------------ //
// CLI parameters
// ------------------------------------------------------------------ //
struct CompareParameters {
    std::vector<std::string> vcfFiles;      // >= 2 phased VCFs (index 0 = truth)
    std::vector<std::string> datasetNames;  // one name per VCF
    std::string              sample;        // empty => first sample in each VCF
    bool                     ignoreSampleName = false;
    bool                     onlySnvs         = false;
    std::string              outPrefix        = "compare_result";
    int                      numThreads       = 1;
    std::string              swBedFile;       // if non-empty, write switch
                                              // error BED to this path

    std::string              version;
    std::string              command;
};

// ------------------------------------------------------------------ //
// One phased heterozygous variant (diploid)
// ------------------------------------------------------------------ //
struct VariantRecord {
    int32_t     pos   = -1;     // 0-based
    std::string ref;
    std::string alt;            // first ALT allele (we are diploid, het)
    // phase: h0, h1 in {0,1}.  -1 means un-phased.
    int8_t      h0 = -1;
    int8_t      h1 = -1;
    int32_t     ps = 0;         // phase-set id (0 = not phased)
    bool        isSnv = true;   // true iff len(ref)==1 && len(alt)==1
};

// Map from variant key (position, ref, alt) to VariantRecord.
// We key on position+ref+alt so the two VCFs can be intersected.
struct VariantKey {
    int32_t     pos;
    std::string ref;
    std::string alt;
    bool operator<(const VariantKey& o) const {
        if (pos != o.pos) return pos < o.pos;
        if (ref != o.ref) return ref < o.ref;
        return alt < o.alt;
    }
    bool operator==(const VariantKey& o) const {
        return pos == o.pos && ref == o.ref && alt == o.alt;
    }
};

// ------------------------------------------------------------------ //
// One BED record for a detected switch error (useful for debugging)
// ------------------------------------------------------------------ //
struct SwBedRecord {
    std::string chromosome;
    int32_t     start = 0;    // 0-based, inclusive  (= position of variant i)
    int32_t     end   = 0;    // 0-based, exclusive  (= position of variant i+1)
    bool        isSnv = true;
    int32_t     psTruth = 0;
    int32_t     psQuery = 0;

    bool operator<(const SwBedRecord& o) const {
        if (chromosome != o.chromosome) return chromosome < o.chromosome;
        if (start      != o.start)      return start      < o.start;
        return end < o.end;
    }
};

// ------------------------------------------------------------------ //
// Per-chromosome comparison result
// ------------------------------------------------------------------ //
struct ChromResult {
    std::string chromosome;

    // --- common-variant statistics --------------------------------- //
    int  commonHet     = 0;   // # het variants common to both VCFs
    int  phasedSnv     = 0;   // # SNVs phased in BOTH data sets
    int  totalSnvTruth = 0;   // # het SNVs in truth (denom. for Phased_SNV%)
    int  phasedIndel     = 0;
    int  totalIndelTruth = 0;

    // --- switch / Hamming (on common & phased variants only) ------- //
    int  intersectionBlocks = 0;       // non-singleton intersection blocks
    int  coveredVariants    = 0;       // variants in those blocks
    int  assessedPairs      = 0;       // sum(block_len-1)
    int  allSwitches        = 0;       // total switch errors
    int  snvSwitches        = 0;       // switch errors restricted to SNV
    int  snvAssessedPairs   = 0;       // pairs restricted to SNV
    int  blockwiseHamming   = 0;       // sum of per-block min-Hamming
    int  totalCompared      = 0;       // variants in the assessed blocks
    int  switchFlipsSwitches = 0;      // switch-flip decomposition
    int  switchFlipsFlips    = 0;

    // --- block statistics on the QUERY (VCF[1]) -------------------- //
    int       numBlocks = 0;           // non-singleton blocks in query
    long long blockSum  = 0;           // sum of block spans in bp
    long long blockN50  = 0;           // N50 of block spans in bp
    std::vector<long long> blockSpans; // raw block spans, used for a
                                       // strict global N50 in the summary

    // --- per-chromosome switch-error BED records ------------------ //
    std::vector<SwBedRecord> swBedRecords;
};

// ------------------------------------------------------------------ //
// CompareProcess – the driver
// ------------------------------------------------------------------ //
class CompareProcess {
public:
    explicit CompareProcess(const CompareParameters& p) : params(p) {}
    int run();

private:
    CompareParameters params;

    // Read all het phased variants from one VCF, keyed by chromosome.
    // Returns map: chrom -> ordered vector<VariantRecord>.
    // Also returns the total number of phased het SNVs in the VCF
    // (regardless of whether they intersect with another VCF).
    bool loadVcf(const std::string& path,
                 std::map<std::string, std::vector<VariantRecord>>& out,
                 std::string& sampleUsed,
                 long long& outTotalPhasedSnv);

    // Compare a single chromosome between truth (idx 0) and query (idx 1).
    // The block-level switch / hamming work is farmed out to threads.
    ChromResult compareChromosome(
            const std::string& chrom,
            const std::vector<VariantRecord>& truth,
            const std::vector<VariantRecord>& query);

    // Compute block N50 given a vector of block spans (bp).
    static long long computeN50(std::vector<long long> spans);

    // Write the final TSV report.
    void writeReport(const std::vector<ChromResult>& chromResults,
                     const std::string& outPath,
                     const std::string& queryBasename,
                     long long queryTotalPhasedSnv,
                     long long queryTotalBlocks,
                     long long queryTotalBlockN50,
                     long long queryTotalBlockSum);
};

#endif
