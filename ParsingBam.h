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
  void writeLine(std::string &input, bool &ps_def,
                 std::ofstream &resultVcf,
                 PhasingResult &phasingResult);

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
  void writeFields(std::ofstream &resultVcf);
  // Emit the three ##INFO header lines for PE / H1 / H2 exactly once.
  // Idempotent: repeated calls are no-ops after the first successful call.
  void writePeInfoHeaders(std::ofstream &resultVcf);
  virtual void writeMetaHeader(const std::string &input, bool &ps_def,
                               std::ofstream &resultVcf);
  void writeColumnHeader(const std::string &input, bool &ps_def,
                         std::ofstream &resultVcf);
  virtual void writeDataLine(const std::string &input,
                             std::ofstream &resultVcf,
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

  // writeDataLine lookups for the chromosome of the last line
  std::string writeChr;
  bool writeChrSet = false;
  std::map<int, SnpVariant> *writeChrVariants = nullptr;
  std::map<int, SnpVariant>::iterator writeVariantCursor;
  const std::set<int> *writeChrFiltered = nullptr;
  // phased positions of each chromosome, sorted, from the PhasingResult
  std::map<std::string, std::vector<std::pair<int, const PhasingElement *>>> writePhased;
  const std::vector<std::pair<int, const PhasingElement *>> *writeChrPhased = nullptr;
  void selectWriteChr(const std::string &chr, PhasingResult &phasingResult);

  // override input parser
  void parserProcess(std::string &input);
  // override output parser
  void writeMetaHeader(const std::string &input, bool &ps_def,
                       std::ofstream &resultVcf);
  void writeDataLine(const std::string &input, std::ofstream &resultVcf,
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

class SVParser : public BaseVairantParser {

private:
  SnpParser *snpFile;

  // chr , variant position (0-base), read
  std::map<std::string, std::map<int, std::map<int, bool>>> *chrVariant;
  // chr, variant position (0-base)
  std::map<std::string, std::map<int, bool>> posDuplicate;

  // override input parser
  void parserProcess(std::string &input);
  // override output parser
  void writeMetaHeader(const std::string &input, bool &ps_def,
                       std::ofstream &resultVcf);
  void writeDataLine(const std::string &input, std::ofstream &resultVcf,
                     PhasingResult &phasingResult);

public:
  SVParser(PhasingParameters &params, SnpParser &snpFile);
  ~SVParser();

  std::map<int, std::map<int, bool>> getVariants(std::string chrName);

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
                       std::ofstream &resultVcf);
  void writeDataLine(const std::string &input, std::ofstream &resultVcf,
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
  std::map<std::string, std::vector<std::pair<int, int>>> SV_map;
  // mod map and iter (the METHParser's map of this chromosome, read in place)
  std::map<int, std::map<std::string, ModVariant>> *currentMod;
  std::map<int, std::map<std::string, ModVariant>> emptyMod;
  std::map<int, std::map<std::string, ModVariant>>::iterator firstModIter;
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
      std::vector<std::pair<int, int>>::iterator &currentSVIter, int &svPos,
      ReadVariant &tmpReadResult);
  void get_snp(const bam_hdr_t &bamHdr, const bam1_t &aln,
               std::vector<ReadVariant> &readVariantVec, ClipCount &clipCount,
               bool hasReference, bool isONT, int svWindow,
               double svThreshold);
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

  void direct_detect_alleles(int lastSNPPos, htsThreadPool &threadPool,
                             PhasingParameters params,
                             std::vector<ReadVariant> &readVariantVec,
                             ClipCount &clipCount,
                             bool hasReference);
};

#endif
