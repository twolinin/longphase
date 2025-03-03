#ifndef UTIL_H
#define UTIL_H

#include <string>
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
#include <iomanip>

struct PhasingElement{
    // i.e. 0|1  or   1|0
    std::string RAstatus;
    int block;
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


// use for parsing
struct Variant{
    Variant(int position, int allele, int quality):
    position(position), 
    allele(allele), 
    quality(quality){};
    
    int position;
    int allele;
    int quality;
    bool underHomopolymer;
};

struct ReadVariant{
    // init function
    ReadVariant(): read_name(""), 
            mapping_quality(0), 
            source_id(""), 
            sample_id(""), 
            reference_start(0), 
            BX_tag(""){}
            
    std::string read_name;
    int mapping_quality;
    std::string source_id;
    std::string sample_id;
    int reference_start;
    std::string BX_tag;
    bool is_reverse;
    std::map<int, int> cnv_mmrate_map;

    std::vector<Variant> variantVec;
    
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

int homopolymerLength(int snp_pos, const std::string &ref_string);

template <typename T>
double calculateMean(const std::vector<T>& data);

#endif