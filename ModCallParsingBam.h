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
#include <deque>
#include <type_traits>

// genotype of a methylation position; "" until it is judged
enum MethGenotype : uint8_t { METH_GT_NONE, METH_GT_HET, METH_GT_HOM_MOD, METH_GT_HOM_UNMOD };
const char *methGenotypeString(uint8_t gt);

struct MethPosInfo {
  MethPosInfo()
      : methreadcnt(0), noisereadcnt(0), canonreadcnt(0), depth(0),
        heterstatus(METH_GT_NONE), strand(-1), variantType(VariantType::MOD) {}

  int methreadcnt;
  int noisereadcnt;
  int canonreadcnt;
  int depth;
  uint8_t heterstatus; // MethGenotype
  int8_t strand; // 0: + , 1: -
  VariantType variantType;
  // reads in the order they were seen, by their index in MethBamParser's name
  // list: index << 1 for a modified call, index << 1 | 1 for an unmodified one
  std::vector<uint32_t> readCalls;
};

// The methylation positions of a chromosome by position, as a std::map<int,
// MethPosInfo> would hold them, in blocks of 1024 bases that are sorted
// arrays. Reads come in order of their start, so the positions are mostly
// added at the end of a block, and the blocks behind the reads can be
// trimmed to their size (shrinkBefore).
class MethPosTable {
public:
  typedef std::pair<int, MethPosInfo> value_type;

private:
  typedef std::vector<value_type> Block;
  static const int blockShift = 10;
  // blocks[i] holds the positions p with p >> blockShift == firstBlock + i
  std::deque<Block> blocks;
  long firstBlock = 0;
  // blocks before this index have been trimmed
  size_t trimmed = 0;

  static long blockOf(int pos) { return long(pos) >> blockShift; }
  Block *blockFor(int pos);
  const Block *blockFor(int pos) const {
    return const_cast<MethPosTable *>(this)->blockFor(pos);
  }

public:
  template <bool Const> class Iter {
    typedef typename std::conditional<Const, const MethPosTable, MethPosTable>::type Table;
    Table *table;
    size_t block, index;
    void skipEmpty() {
      while (block < table->blocks.size() && index >= table->blocks[block].size()) {
        block++;
        index = 0;
      }
    }

  public:
    typedef typename std::conditional<Const, const value_type, value_type>::type Value;
    Iter(Table *t, size_t b, size_t i) : table(t), block(b), index(i) { skipEmpty(); }
    Value &operator*() const { return table->blocks[block][index]; }
    Value *operator->() const { return &table->blocks[block][index]; }
    Iter &operator++() {
      index++;
      skipEmpty();
      return *this;
    }
    Iter operator++(int) {
      Iter old = *this;
      ++*this;
      return old;
    }
    bool operator==(const Iter &o) const { return block == o.block && index == o.index; }
    bool operator!=(const Iter &o) const { return !(*this == o); }
  };
  typedef Iter<false> iterator;
  typedef Iter<true> const_iterator;

  iterator begin() { return iterator(this, 0, 0); }
  iterator end() { return iterator(this, blocks.size(), 0); }
  const_iterator begin() const { return const_iterator(this, 0, 0); }
  const_iterator end() const { return const_iterator(this, blocks.size(), 0); }
  iterator find(int pos);
  const_iterator find(int pos) const;
  // the entry at pos, added if absent (as std::map::operator[])
  MethPosInfo &operator[](int pos);
  // trims the blocks that lie wholly before pos
  void shrinkBefore(int pos);
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

  MethPosTable *chrMethMap;
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
  // every variant position, include SNP and SV, sorted by position: the
  // number of distinct read names that cover it, and the type recorded by the
  // read whose name sorts first (the last record of that name wins)
  struct NodeSummary {
    int readCount;
    int type;
  };
  typedef std::vector<std::pair<int, NodeSummary>> NodeList;
  NodeList nodeInfo;
  NodeList::const_iterator findNode(int position) const;
  int nodeReadCount(int position) const;
  // check the variant type of the position
  int checkVariantType(int position);
  float computeMinConnection(int pos1, int pos2) const;
  void buildNodeInfo(const std::vector<ModReadVariant> &reads);
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
