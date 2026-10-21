#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <vector>
#include <omp.h>
#include <numeric> 
#include <set>
#include <functional>
#include <htslib/vcf.h>

struct PhasingElement{
    // i.e. 0|1  or   1|0
    std::string RAstatus;
    int block;
    // vote-based diagnostics attached only when the position was phased.
    // Consumers (writeDataLine) check whether the element exists in
    // PhasingResult; if it does, these three values are always valid.
    float h1 = 0.0f;      // weighted HP1 vote count at this position
    float h2 = 0.0f;      // weighted HP2 vote count at this position
    float entropy = 0.0f; // Shannon entropy of (h1, h2), in bits (0..1)
};
typedef std::map<std::string,PhasingElement> PhasingResult;
typedef std::map<std::string,PhasingResult> ChrPhasingResult;

/**
 * Merges phasing results from all chromosomes into one result.
 *
 * @param allChrPhasingResults Contains phasing results for each chromosome.
 * @param mergedPhasingResult A map to store the combined phasing results.
 *
 * The function combines the phasing elements from each chromosome into mergedPhasingResult.
 */
void mergeAllChrPhasingResult(const ChrPhasingResult& allChrPhasingResults, PhasingResult& mergedPhasingResult);

/**
 * Configures thread counts for chromosome processing and BAM parsing.
 *
 * @param defaultThreads Initial thread count for chromosome processing.
 * @param availableThreads  Total threads available for processing
 * @param chrnumThreads Set to the number of threads for chromosome processing.
 * @param bamParsernumThreads Set to the number of threads for BAM parsing.
 */
void setPhasingNumThreads(const int& defaultChrThreads,const int& availableThreads,  int& chrnumThreads, int& bamParsernumThreads);

/**
 * Configures thread counts for chromosome processing and BAM parsing.
 *
 * @param availableThreads  Total threads available for processing
 * @param chrnumThreads Set to the number of threads for chromosome processing.
 * @param bamParsernumThreads Set to the number of threads for BAM parsing.
 */
void setModcallNumThreads(const int& availableThreads,  int& chrNumThreads, int& bamParserNumThreads);

enum VariantType : int8_t {
    MOD = 0,
    SNP = 1,
    INDEL = 2,
    SV = 3
};

// use for parsing
struct Variant{
    Variant(int position, int allele, int quality, VariantType type):
    position(position), 
    quality(quality),
    allele(allele), 
    type(type){};
    
    // 8 bytes: reads hold many of these
    int position;
    int16_t quality;     // base quality (0-255), or a negative code for SV / modification / indel
    int8_t allele;
    VariantType type;
};

struct ReadVariant{
    // init function
    ReadVariant(): read_name(""), 
            reference_start(0){}
            
    // 64 bytes: phase holds one per alignment with variants
    std::string read_name;
    std::vector<Variant> variantVec;
    int reference_start;
    bool is_reverse;
    
    void sort();
};

struct less_than_key
{
    inline bool operator() (const Variant& v1, const Variant& v2)
    {
        return (v1.position < v2.position);
    }
};


std::string getTargetString(std::string line, std::string start_sign, std::string end_sign);

// Returns the character start position in `sample` for the value of `tag` in `format`.
int findFieldValueStart(const std::string &format, const std::string &sample,
                        const std::string &tag);

// Marks indels in tandem repeats (2-mer repeated >= 5 times) as dangerous.
bool isDangerIndel(int variant_pos, const std::string &ref, size_t ref_allele_len,
                   size_t alt_allele_len);

int homopolymerLength(int snp_pos, const std::string &ref_string);

template <typename T>
double calculateMean(const std::vector<T>& data);

// readVcfRecords reads the records of fp as a bcf_read loop does: in file
// order, up to the first one that does not parse. They are read in batches;
// each batch is parsed on up to `threads` threads, each with its own copy of
// the header (parsing adds contigs and tags missing from the header to it),
// and work(i, hdr, rec) is called there for the i-th record of the batch.
// Then done(i, hdr, rec) is called for the records in file order. work may
// only write data of record i.
typedef std::function<void(size_t, bcf_hdr_t*, bcf1_t*)> VcfRecordFn;
void readVcfRecords(htsFile* fp, bcf_hdr_t* hdr, int threads, size_t batch,
                    const VcfRecordFn& work, const VcfRecordFn& done);

#endif
