#ifndef MODCALLPARSINGBAM_H
#define MODCALLPARSINGBAM_H

#include "ModCallProcess.h"
#include "ParsingBam.h"
#include "PhasingGraph.h"
#include "PhasingProcess.h"
#include "Util.h"
#include <htslib/faidx.h>
#include <htslib/kbitset.h>
#include <htslib/khash.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

struct MethPosInfo {
  MethPosInfo()
      : methreadcnt(0), noisereadcnt(0), canonreadcnt(0), depth(0),
        heterstatus(""), strand(-1), variantType(VariantType::MOD) {}

  int methreadcnt;
  int noisereadcnt;
  int canonreadcnt;
  int depth;
  std::string heterstatus;
  int strand; // 0: + , 1: -
  VariantType variantType;
  // reads by their index in MethBamParser's name list
  std::vector<uint32_t> modReadVec;
  std::vector<uint32_t> nonModReadVec;
};

// The ReferenceChromosome struct is used to store information about a reference
// fasta.
struct ReferenceChromosome {
  std::string name;
  int length;

  ReferenceChromosome(const std::string &name, int length)
      : name(name), length(length) {}
};

// FastaParser: to get each chromosome name and lastpos(for chunksize)
class MethFastaParser {
private:
  // file name
  std::string fastaFile;

public:
  MethFastaParser(std::string fastaFile,
                  std::vector<ReferenceChromosome> &chrInfo);
  ~MethFastaParser();
  static std::string fetchSequence(const std::string &fastaFile,
                                   const ReferenceChromosome &chr);
};

// A chromosome's reference in about a quarter of the memory: A, C, G and T in
// 2 bits each, with the lowercase stretches and the runs of any other
// character listed apart. at() gives back the original character.
class PackedReference {
public:
  explicit PackedReference(const std::string &seq);
  size_t length() const { return len; }
  char at(size_t pos) const;

private:
  struct Run {
    uint32_t start, end; // [start, end)
    char c;
  };
  size_t len;
  std::vector<uint8_t> bases;  // 2 bits per position, 0 where the character is in `other`
  std::vector<Run> lowercase;  // runs of a, c, g, t (c unused)
  std::vector<Run> other;      // runs of one character that is not A, C, G, T, a, c, g, t
  static const Run *findRun(const std::vector<Run> &runs, size_t pos);
};

struct MethPosProb {
  MethPosProb(int position, int prob) : position(position), prob(prob) {};

  int position;
  int prob;
};

struct AlignmentMethExtend : Alignment {
  // hts_base_mod_state *m;
  // store the read position and probability
  std::vector<MethPosProb> queryMethVec;
};

class MethSnpParser : public BaseVairantParser {
private:
  ModCallParameters *params; // Use ModCallParameters instead
  // std::map<std::string, std::map<int, SnpVariant>> *chrVariant;
  bool hasSnpData; // indicate if the SNP data is loaded successfully

  bool loadSNPVCF();
  void parserProcess(std::string &input);
  void writeDataLine(const std::string &input, std::ostream &resultVcf,
                     PhasingResult &phasingResult);

public:
  MethSnpParser(ModCallParameters &in_params); // New constructor
  ~MethSnpParser();

  std::map<std::string, std::map<int, SnpVariant>> *chrVariant;

  std::map<int, SnpVariant> getVariants_markindel(std::string chrName,
                                              const std::string &ref);
  void writeResult(PhasingResult phasingResult);
  bool hasValidSnpData() const; // check if the SNP data is loaded successfully
};

class MethBamParser {
private:
  ModCallParameters *params;
  std::string chrName;
  int refstartpos;
  // what the reads and the output need of the chromosome's reference
  PackedReference reference;

  std::map<int, MethPosInfo> *chrMethMap;
  // names of the reads parsed so far, one after another; read i starts at nameStart[i]
  std::string readNames;
  std::vector<size_t> nameStart;
  //<pos, <positive strand cnt, negative strand cnt>>
  std::map<int, std::pair<int, int>> *readStartEndMap;

  // get methylation tag from BAM file
  std::string get_aux_tag(const bam1_t *aln, const char tag[2]);
  bool isValidAlignment(const bam1_t *aln) const;
  bool extractSnpAllele(
      std::map<int, SnpVariant>::iterator &currentVariantIter, int ref_pos,
      int length, int querypos, int cigaridx, int aln_core_n_cigar,
      const uint32_t *cigar, const bam1_t &aln, ModReadVariant *tmpReadResult);
  void classifyMethylationBase(int &n, int &pos, hts_base_mod *mods,
                               int querypos, int length, int refpos,
                               const bam1_t &aln, ModReadVariant *tmpReadResult,
                               hts_base_mod_state *mod_state);
  bool handleDeletionInHomopolymer(
      std::map<int, SnpVariant>::iterator &currentVariantIter, int ref_pos,
      int length, int querypos, const bam1_t &aln,
      hts_base_mod_state *mod_state, ModReadVariant *tmpReadResult);
  void updateDepthBoundary(int refstart, int refpos, bool is_reverse);
  void computeSinglePositionGenotype();
  std::set<int> applyCpGStrandPairGenotype();
  void filterModReadVariants(const std::set<int> &positionPairs,
                             std::vector<ModReadVariant> &readVariantVec,
                             std::vector<ModReadVariant> &modReadVariantVec);
  std::string buildInfoStr(const MethPosInfo &info) const;
  std::string formatMethVcfLine(int pos, const MethPosInfo &info,
                                int chrLen) const;
  void parse_CIGAR(const bam_hdr_t &bamHdr, const bam1_t &aln,
                   std::vector<ModReadVariant> &readVariantVec);

  // record methylation position and probability from methylation tag
  void getmeth(AlignmentMethExtend &align);
  std::map<int, SnpVariant> *currentVariants;
  std::map<int, SnpVariant>::iterator firstVariantIter;

public:
  MethBamParser(std::string inputChrName, ModCallParameters &in_params,
                MethSnpParser &snpMap, const std::string &ref_string);
  ~MethBamParser();
  void detectMeth(std::string chrName, int lastSNPPos,
                  htsThreadPool &threadPool,
                  std::vector<ModReadVariant> &readVariantVec);
  void dumpMethMap(const std::string &chrName, std::ofstream &out) const;
  void dumpReadStartEnd(const std::string &chrName, std::ofstream &out) const;
  void calculateDepth();
  void exportResult(const std::string &chrName, int chrLen,
                    std::vector<int> &passPosition,
                    std::ostringstream &modCallResult);
  void judgeMethGenotype(std::string chrName,
                         std::vector<ModReadVariant> &readVariantVec,
                         std::vector<ModReadVariant> &modReadVariantVec);
  void dumpMethGenoMap(const std::string &chrName, std::ofstream &out) const;
};

void writeResultVCF(ModCallParameters &params,
                    std::vector<ReferenceChromosome> &chrInfo,
                    std::map<std::string, std::ostringstream> &chrResult);

class MethylationGraph {
private:
  ModCallParameters *params;
  std::vector<ModReadVariant> *readVariant;

  // modifications
  // position, quality
  std::map<int, ReadBaseMap *> *forwardModNode;
  std::map<int, ReadBaseMap *> *reverseModNode;

  // By default, a Map in C++ is sorted in increasing order based on its key.
  // position, edge
  std::map<int, VariantEdge *> *edgeList;
  // record all variant position, include SNP and SV
  // <position, <readID, VariantType>>
  std::map<int, std::map<std::string, VariantType>> *nodeInfo;
  // check the variant type of the position
  int checkVariantType(int position);
  float computeMinConnection(int pos1, int pos2) const;
  double computeMajorRatio(std::pair<int, int> counts) const;
  bool isGoodConnection(std::pair<int, int> counts, float minConn) const;
  void collectAllModPoints(std::set<int> &strongPoints);
  void collectStrongPoints(std::set<int> &strongPoints,
                           std::set<int> &weakPoints,
                           std::vector<int> &prepass);
  void expandStrongPoints(const std::set<int> &strongPoints,
                          std::vector<int> &prepass,
                          std::set<int> &weakPoints,
                          std::set<int> &addedPositions,
                          bool hasValidSnpData);
  void iterateWeakPoints(std::set<int> &weakPoints,
                         std::set<int> &addedPositions,
                         std::vector<int> &prepass);
  void filterByNeighborConnection(std::vector<int> &prepass,
                                  std::vector<int> &passPosition);

public:
  MethylationGraph(ModCallParameters &params);
  ~MethylationGraph();

  void addEdge(std::vector<ModReadVariant> &in_readVariant, std::string chrName);
  void connectResults(std::string chrName, std::vector<int> &passPosition,
                      bool hasValidSnpData = true);
  void destroy();
};

#endif
