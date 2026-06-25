#ifndef HAPLOTAGPROCESS_H
#define HAPLOTAGPROCESS_H

#include "Util.h"
#include "ParsingBam.h"
#include <cmath>

struct HaplotagParameters
{
    int numThreads;
    int qualityThreshold;
    int svWindow;
    
    double percentageThreshold;
    double svThreshold;
    
    std::string snpFile;
    std::string svFile;
    std::string modFile;
    std::string bamFile;
    std::string fastaFile;
    std::string resultPrefix;
    std::string region;
    std::string command;
    std::string version;
    std::string outputFormat;
    
    bool tagSupplementary;
    bool writeReadLog;
};

class HaplotagProcess
{
    void variantParser(std::string variantFile);
    void compressParser(std::string &variantFile);
    void unCompressParser(std::string &variantFile);
    void parserProcess(std::string &input);
    void parseHeaderLine(const std::string &line);
    void parseRecordLine(const std::string &line);
    void processSnpRecord(const std::vector<std::string> &fields, int pos, int haploDir, const std::string &psValue);
    void processSvRecord(const std::vector<std::string> &fields, int haploDir);
    void processModRecord(const std::vector<std::string> &fields, int haploDir);
    
    std::vector<int> collectLastVariantPos() const;
    std::ofstream* initLogFile(const HaplotagParameters &params);
    void resolveRegionScope(const HaplotagParameters &params);
    void tagChromosome(const std::string &chr, const HaplotagParameters &params,
                       samFile *in, samFile *out, bam_hdr_t *bamHdr, hts_idx_t *idx,
                       bam1_t *aln, const std::string &chr_reference, std::ofstream *tagResult);
    void writeUnmappedTail(samFile *in, samFile *out, bam_hdr_t *bamHdr, hts_idx_t *idx);
    void tagRead(HaplotagParameters &params);
    
    std::vector<std::string> chrVec;
    std::map<std::string, int> chrLength;
    
    // chr, variant position (0-base), allele haplotype
    std::map<std::string, std::map<int, SnpVariant> > chrVariant;
    // chr, variant position (0-base), phased set
    std::map<std::string, int > psIndex;
    std::map<std::string, std::map<int, int> > chrVariantPS;
    
    std::map<std::string, std::map<int, std::string> > chrVariantHP1;
    std::map<std::string, std::map<int, std::string> > chrVariantHP2;

    std::map<int, SnpVariant> currentChrVariants;
    std::map<int, SnpVariant>::iterator firstVariantIter;
    // The number of SVs occurring on different haplotypes in a read
    std::map<std::string, std::map<int, int> > readSVHapCount;

    // record SV region
    std::map<std::string, std::vector<std::tuple<int, int, int> > > chrRegions;
    std::vector<std::tuple<int, int, int> > currentchrRegions;
    std::vector<std::tuple<int, int, int> >::iterator firstSVIter;

    void initFlag(bam1_t *aln, std::string flag);

    void collectSVWindowEvidence(const bam1_t &aln, int cigarIdx, int n_cigar,
                                 const uint32_t *cigar,
                                 std::vector<std::tuple<int, int, int> >::iterator svIter,
                                 int svWindow, double svThreshold);
    void collectMatchOpEvidence(const bam1_t &aln, const std::string &chrName,
                                int ref_pos, int query_pos, int length,
                                int cigarIdx, int n_cigar, const uint32_t *cigar,
                                std::map<int, SnpVariant>::iterator &variantIter,
                                std::map<int, int> &variantsHP, std::map<int, int> &countPS,
                                int &hp1Count, int &hp2Count);
    void collectDeletionOpEvidence(const bam1_t &aln, const std::string &chrName,
                                   int ref_pos, int query_pos, int length,
                                   const std::string &ref_string,
                                   std::map<int, SnpVariant>::iterator &variantIter,
                                   std::map<int, int> &variantsHP, std::map<int, int> &countPS,
                                   int &hp1Count, int &hp2Count);
    int decideHaplotype(int hp1Count, int hp2Count, const std::map<int, int> &countPS,
                        double percentageThreshold, int &pqValue);
    void writeReadLogRow(const bam_hdr_t &bamHdr, const bam1_t &aln,
                         int hpResult, int hp1Count, int hp2Count,
                         const std::map<int, int> &variantsHP, const std::map<int, int> &countPS,
                         int pqValue, std::ofstream *tagResult);
    int judgeHaplotype(const  bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, const std::string &ref_string, int svWindow, double svThreshold);
    
    int totalAlignment;
    int totalSupplementary;
    int totalSecondary;
    int totalUnmapped;
    int totalTagCount;
    int totalUnTagCount;
    
    std::time_t processBegin;
    bool integerPS;
    bool parseSnpFile;
    bool parseSVFile;
    bool parseMODFile;
    public:
        HaplotagProcess(HaplotagParameters params);
        ~HaplotagProcess();

};


#endif
