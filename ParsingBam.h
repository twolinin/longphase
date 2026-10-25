#ifndef PARSINGBAM_H
#define PARSINGBAM_H

#include "ModCallProcess.h"
#include "PhasingProcess.h"
#include "Util.h"
#include <htslib/faidx.h>
#include <htslib/kbitset.h>
#include <htslib/khash.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <functional>
#include <memory>
#include <unordered_set>
#include <zlib.h>

struct SnpVariant {
  std::string Ref;
  std::string Alt;
  bool is_danger = false;
  // homopolymerLength() at this position, set with is_danger while the
  // chromosome's reference is loaded (markChrVariants)
  int8_t homopolymer = 1;
};

struct ModVariant {
  bool is_reverse = false;
  bool is_modify = false;
};


// Reads sequences of a FASTA file as faidx_fetch_seq does. htslib reads the
// file one character at a time; for an uncompressed FASTA with its .fai and
// regular lines this reads the same bytes in blocks. Otherwise (compressed
// FASTA, no .fai, irregular lines) it calls htslib. One reader per thread.
class FastaReader {
public:
  explicit FastaReader(const std::string &fastaFile);
  ~FastaReader();
  bool loaded() const { return fai != nullptr; }
  // faidx_fetch_seq(fai, chr, beg, end, &len) as a string ("" when it returns
  // NULL). With cache, a block of sequence around the range is kept, for many
  // small fetches in increasing order.
  std::string fetch(const std::string &chr, int beg, int end, int &len,
                    bool cache);

private:
  faidx_t *fai;
  int fd = -1; // the FASTA, if it can be read in blocks
  struct Entry {
    int64_t len, offset;
    int lineBases, lineWidth;
    bool irregular = false;
  };
  std::map<std::string, Entry> index;
  // cached bytes of the file, from rawStart
  int64_t rawStart = -1;
  std::vector<char> raw;
  // sequence [beg, end] of e into out, reading the file in blocks or from
  // the cached bytes; false if the lines are not as the .fai describes or
  // the file is short
  bool readSeq(const Entry &e, int64_t beg, int64_t end, std::string &out,
               bool cache);
  FastaReader(const FastaReader &) = delete;
  FastaReader &operator=(const FastaReader &) = delete;
};

// The reference of chr up to 5 bases past lastPos (0-based), as
// FastaParser reads it. loadChrReference opens the FASTA itself, so
// several threads can call it at once.
std::string fetchChrReference(faidx_t *fai, const std::string &chr, int lastPos);
std::string loadChrReference(const std::string &fastaFile,
                             const std::string &chr, int lastPos);

class FastaParser {
private:
  // file name
  std::string fastaFile;
  std::vector<std::string> chrName;
  std::vector<int> last_pos;

public:
  FastaParser(std::string fastaFile, std::vector<std::string> chrName,
              std::vector<int> last_pos, int numThreads);
  ~FastaParser();

  // chrName, chr string
  std::map<std::string, std::string> chrString;
};

class BaseVairantParser {
public:
  BaseVairantParser() : params(nullptr), commandLine(false), pe_def(false) {}
  // input parser
  void compressParser(std::string &variantFile);
  void unCompressParser(std::string &variantFile);
  virtual void parserProcess(std::string &input) = 0;
  // output parser
  void compressInput(std::string variantFile, std::string resultFile,
                     PhasingResult &phasingResult);
  void unCompressInput(std::string variantFile, std::string resultFile,
                       PhasingResult &phasingResult);
  void dispatchWriteResult(const std::string &inputFile,
                           const std::string &outputFile,
                           PhasingResult &phasingResult);
  virtual void writeLine(std::string &input, bool &ps_def,
                         std::ostream &resultVcf,
                         PhasingResult &phasingResult);
  // called after the last line
  virtual void finishWriting(std::ostream &resultVcf,
                             PhasingResult &phasingResult) {}
  virtual ~BaseVairantParser() {}

protected:
  PhasingParameters *params;
  bool commandLine;
  // Tracks whether ##INFO=<ID=PE|H1|H2,...> headers have been emitted
  // so the three lines are written exactly once per output VCF, even
  // when the input already declares them.
  bool pe_def;
  // writeDataLine buffers, reused from line to line
  std::vector<std::string> lineFields;
  std::string lineOut;
  // join fields with tabs and write them as one line
  void writeFields(std::ostream &resultVcf);
  // Emit the three ##INFO header lines for PE / H1 / H2 exactly once.
  // Idempotent: repeated calls are no-ops after the first successful call.
  void writePeInfoHeaders(std::ostream &resultVcf);
  virtual void writeMetaHeader(const std::string &input, bool &ps_def,
                               std::ostream &resultVcf);
  void writeColumnHeader(const std::string &input, bool &ps_def,
                         std::ostream &resultVcf);
  virtual void writeDataLine(const std::string &input,
                             std::ostream &resultVcf,
                             PhasingResult &phasingResult) = 0;

private:
  void readGzLines(const std::string &variantFile,
                   std::function<void(const std::string &)> callback);
};

class SnpParser : public BaseVairantParser {

private:
  // std::string variantFile;
  //  chr, variant position (0-base), allele haplotype
  std::map<std::string, std::map<int, SnpVariant>> *chrVariant;
  // id and idx
  std::vector<std::string> chrName;
  // chr, variant position (0-base)
  std::map<std::string, std::map<int, bool>> chrVariantHomopolymer;

  // Track the position of filtered indels
  std::map<std::string, std::set<int>> filteredIndelPositions;
  // <prefix>_removed_indels.log, open only when --indelQuality is set
  std::ofstream removedIndelsLog;
  // column of the phased sample (--sample) in a VCF data line
  size_t sampleCol = 9;

  // What writing a data line looks up for its chromosome: its SNPs, its
  // filtered indels and its phased positions (sorted, from the PhasingResult)
  struct WriteChr {
    const std::map<int, SnpVariant> *variants;
    const std::set<int> *filtered;
    std::vector<std::pair<int, const PhasingElement *>> phased;
  };
  std::map<std::string, WriteChr> writeChrs;
  const WriteChr &writeChrLookup(const std::string &chr,
                                 PhasingResult &phasingResult);
  // The data lines are formatted in batches on several threads and written
  // in order; a header line first writes the data lines before it.
  std::vector<std::string> pendingLines;
  size_t pendingCount = 0;
  std::vector<const WriteChr *> pendingChr;
  std::vector<std::string> pendingOut;
  std::vector<std::vector<std::string>> threadFields;
  std::unique_ptr<BatchThreads> writeThreads;
  void formatDataLine(const std::string &input, const WriteChr &chr,
                      std::vector<std::string> &fields,
                      std::string &out) const;
  void flushDataLines(std::ostream &resultVcf, PhasingResult &phasingResult);
  void writeLine(std::string &input, bool &ps_def, std::ostream &resultVcf,
                 PhasingResult &phasingResult) override;
  void finishWriting(std::ostream &resultVcf,
                     PhasingResult &phasingResult) override;

  // override input parser
  void parserProcess(std::string &input);
  // override output parser
  void writeMetaHeader(const std::string &input, bool &ps_def,
                       std::ostream &resultVcf);
  void writeDataLine(const std::string &input, std::ostream &resultVcf,
                     PhasingResult &phasingResult);

public:
  SnpParser(PhasingParameters &in_params);
  ~SnpParser();

  std::map<int, SnpVariant> getVariants(std::string chrName);

  // Set is_danger and homopolymer of chrName's variants from its reference.
  // These are the only uses of the reference while phasing, so it can be
  // freed before the reads are scanned.
  void markChrVariants(const std::string &chrName, const std::string &ref);

  // chrName's variants, or nullptr
  std::map<int, SnpVariant> *chrVariants(const std::string &chrName);

  std::vector<std::string> getChrVec();

  int getLastSNP(std::string chrName);

  void writeResult(PhasingResult &phasingResult, const std::string &outPrefix);

  bool findSNP(std::string chr, int posistion);

  void filterSNP(std::string chr, std::vector<ReadVariant> &readVariantVec);
};

// Read support of an SV taken from the SV VCF, used next to the CIGAR to decide
// which reads carry the SV (issue #126).
struct SVReadSupport {
  // from SVTYPE, or from the sign of SVLEN when SVTYPE is absent
  bool isDeletion = false;
  // INFO/RNAMES (e.g. Sniffles2 --output-rnames): reads supporting the SV
  bool hasReadNames = false;
  std::unordered_set<std::string> readNames;
};

// For each deletion of the SV VCF, how far its start can move in a repeat of
// the reference without changing the sequence left by the deletion (issue
// #126 compares a read's split junction with it). Worked out while the
// chromosome's reference is held, since phase frees it before the reads.
class SVDeletionRepeats {
public:
  // the deletions among svs (position, SVLEN; as SVParser::getVariants gives
  // them), anchored at their position and at the base before it
  void build(const std::map<int, std::map<int, bool>> &svs,
             const std::string &ref);
  // as deleting ref[a, a + length) and ref[b, b + length) leaving the same
  // sequence, for b and length given to build (b the position or the base
  // before it) and |a - b| <= max(maxShift, length) + 1
  bool sameDeletion(int a, int b, int length) const;
  // the largest shift of a deletion start that is asked about
  static const int maxShift = 1000;

private:
  struct Anchor {
    int pos;
    int length;
    // the deletion at any start in [low, high] leaves the same sequence
    int low;
    int high;
    bool operator<(const Anchor &o) const {
      return pos != o.pos ? pos < o.pos : length < o.length;
    }
  };
  std::vector<Anchor> anchors;
  size_t refLength = 0;
};

class SVParser : public BaseVairantParser {

private:
  SnpParser *snpFile;

  // chr , variant position (0-base), read
  std::map<std::string, std::map<int, std::map<int, bool>>> *chrVariant;
  // chr, variant position (0-base)
  std::map<std::string, std::map<int, bool>> posDuplicate;
  // chr, variant position (1-base, same key as chrVariant), read support
  std::map<std::string, std::map<int, SVReadSupport>> chrSVSupport;
  // chr, (position (1-base), |SVLEN|) of every deletion record of the SV VCF,
  // whatever its genotype and whether phase uses it, sorted
  std::map<std::string, std::vector<std::pair<int, int>>> chrDeletionRecords;

  // override input parser
  void parserProcess(std::string &input);
  // override output parser
  void writeMetaHeader(const std::string &input, bool &ps_def,
                       std::ostream &resultVcf);
  void writeDataLine(const std::string &input, std::ostream &resultVcf,
                     PhasingResult &phasingResult);

public:
  SVParser(PhasingParameters &params, SnpParser &snpFile);
  ~SVParser();

  std::map<int, std::map<int, bool>> getVariants(std::string chrName);
  // read support of the SVs on chrName; nullptr if there is none
  const std::map<int, SVReadSupport> *
  getReadSupport(const std::string &chrName) const;
  // deletion records on chrName (see chrDeletionRecords); nullptr if none
  const std::vector<std::pair<int, int>> *
  getDeletionRecords(const std::string &chrName) const;

  void writeResult(PhasingResult &phasingResult, const std::string &outPrefix);

  bool findSV(std::string chr, int posistion);
};

class METHParser : public BaseVairantParser {

private:
  SnpParser *snpFile;
  SVParser *svFile;

  int representativePos;
  int upMethPos;

  // chr , variant position (0-base), read, (is reverse strand)
  std::map<std::string, std::map<int, std::map<std::string, ModVariant>>>
      *chrVariant;
  // In a series of consecutive methylation positions,
  // the first methylation position will be used as the representative after
  // merging. This map is used to find the coordinates that originally
  // represented itself
  std::map<int, int> *representativeMap;

  // override input parser
  void parserProcess(std::string &input);
  // override output parser
  void writeMetaHeader(const std::string &input, bool &ps_def,
                       std::ostream &resultVcf);
  void writeDataLine(const std::string &input, std::ostream &resultVcf,
                     PhasingResult &phasingResult);

public:
  std::map<int, std::map<std::string, ModVariant>>
  getVariants(std::string chrName);

  // chrName's modifications, or nullptr
  std::map<int, std::map<std::string, ModVariant>> *
  chrVariants(const std::string &chrName);

  METHParser(PhasingParameters &params, SnpParser &snpFile, SVParser &svFile);
  ~METHParser();

  void writeResult(PhasingResult &phasingResult, const std::string &outPrefix);
};

struct Alignment {
  std::string chr;
  std::string qname;
  int refStart;
  int qlen;
  char *qseq;
  int cigar_len;
  int *op;
  int *ol;
  char *quality;
  bool is_reverse;
};

enum ClipFrontBack { FRONT = 0, BACK = 1 };

// pos<read start|read end ,count >
using ClipCount = std::map<int, std::map<ClipFrontBack, int>>;

class BamParser {

private:
  std::string chrName;
  std::vector<std::string> BamFileVec;
  // SNP map and iter (the SnpParser's map of this chromosome, read in place)
  std::map<int, SnpVariant> *currentVariants;
  std::map<int, SnpVariant>::iterator firstVariantIter;
  // SV map and iter
  std::map<int, std::map<int, bool>> *currentSV;
  std::vector<std::pair<int, int>>::iterator firstSVIter;
  // the chromosome's SVs: (start, length)
  std::vector<std::pair<int, int>> svList;
  // read support of the SVs on this chromosome (owned by SVParser)
  const std::map<int, SVReadSupport> *currentSVSupport = nullptr;
  // deletion records on this chromosome (owned by SVParser)
  const std::vector<std::pair<int, int>> *currentDeletionRecords = nullptr;
  // repeats around the chromosome's SV deletions (made from the reference)
  SVDeletionRepeats deletionRepeats;
  // mod map and iter (the METHParser's map of this chromosome, read in place)
  std::map<int, std::map<std::string, ModVariant>> *currentMod;
  std::map<int, std::map<std::string, ModVariant>> emptyMod;
  std::map<int, std::map<std::string, ModVariant>>::iterator firstModIter;
  // names of the reads kept in readVariantVec, one after another (read_id
  // indexes nameStart until rankReadNames turns it into the name's rank)
  std::string readNames;
  std::vector<size_t> nameStart;
  void rankReadNames(std::vector<ReadVariant> &readVariantVec);
  void initSnpVariants(SnpParser &snpMap);
  void initSvVariants(SVParser &svFile);
  void initModVariants(METHParser &modFile);
  bool handleDeletionHomopolymerSNP(
      const bam1_t &aln, int ref_pos, int query_pos, int del_len,
      std::map<int, SnpVariant>::iterator &currentVariantIter,
      ReadVariant &tmpReadResult);
  void handleSNPVariant(const bam1_t &aln, int cigar_idx, int cigar_oplen,
                        int ref_pos, int query_pos,
                        std::map<int, SnpVariant>::iterator &currentVariantIter,
                        int &variantPos, ReadVariant &tmpReadResult);
  void handleModVariant(
      const bam1_t &aln, int variantPos,
      std::map<int, std::map<std::string, ModVariant>>::iterator &currentModIter,
      int &modPos, ReadVariant &tmpReadResult);
  void handleSVVariant(
      const bam1_t &aln, int cigar_i, int svWindow, double svThreshold,
      int mappingQuality,
      std::vector<std::pair<int, int>>::iterator &currentSVIter, int &svPos,
      ReadVariant &tmpReadResult);
  void get_snp(const bam_hdr_t &bamHdr, const bam1_t &aln,
               std::vector<ReadVariant> &readVariantVec, ClipCount &clipCount,
               bool hasReference, bool isONT, int svWindow,
               double svThreshold, int mappingQuality);
  void getClip(int pos, int clipFrontBack, int len, ClipCount &clipCount);
  bool isValidAlignment(const bam1_t &aln, int mappingQuality);
  void processBamFile(const std::string &bamFile, int lastSNPPos,
                      htsThreadPool &threadPool, const PhasingParameters &params,
                      std::vector<ReadVariant> &readVariantVec,
                      ClipCount &clipCount, bool hasReference);

public:
  BamParser(std::string chrName, std::vector<std::string> inputBamFileVec,
            SnpParser &snpMap, SVParser &svFile, METHParser &modFile);
  ~BamParser();
  // repeats around the chromosome's SV deletions, made while the reference
  // was held (see SVDeletionRepeats)
  void setDeletionRepeats(SVDeletionRepeats repeats) {
    deletionRepeats = std::move(repeats);
  }

  void direct_detect_alleles(int lastSNPPos, htsThreadPool &threadPool,
                             PhasingParameters params,
                             std::vector<ReadVariant> &readVariantVec,
                             ClipCount &clipCount,
                             bool hasReference);
};

#endif
