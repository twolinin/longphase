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
};

struct ModVariant {
  bool is_reverse = false;
  bool is_modify = false;
};


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
  BaseVairantParser() : params(nullptr), commandLine(false) {}
  // input parser
  void compressParser(std::string &variantFile);
  void unCompressParser(std::string &variantFile);
  virtual void parserProcess(std::string &input) = 0;
  // output parser
  void compressInput(std::string variantFile, std::string resultFile,
                     PhasingResult phasingResult);
  void unCompressInput(std::string variantFile, std::string resultFile,
                       PhasingResult phasingResult);
  void dispatchWriteResult(const std::string &inputFile,
                           const std::string &outputFile,
                           PhasingResult phasingResult);
  void writeLine(std::string &input, bool &ps_def,
                 std::ofstream &resultVcf,
                 PhasingResult &phasingResult);

protected:
  PhasingParameters *params;
  bool commandLine;
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

  std::map<int, SnpVariant> getVariants_markindel(std::string chrName,
                                                  const std::string &ref);

  std::vector<std::string> getChrVec();

  int getLastSNP(std::string chrName);

  void writeResult(PhasingResult phasingResult);

  bool findSNP(std::string chr, int posistion);

  void filterSNP(std::string chr, std::vector<ReadVariant> &readVariantVec,
                 std::string &chr_reference);
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
  void writeDataLine(const std::string &input, std::ofstream &resultVcf,
                     PhasingResult &phasingResult);

public:
  SVParser(PhasingParameters &params, SnpParser &snpFile);
  ~SVParser();

  std::map<int, std::map<int, bool>> getVariants(std::string chrName);

  void writeResult(PhasingResult phasingResult);

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
  void writeDataLine(const std::string &input, std::ofstream &resultVcf,
                     PhasingResult &phasingResult);

public:
  std::map<int, std::map<std::string, ModVariant>>
  getVariants(std::string chrName);

  METHParser(PhasingParameters &params, SnpParser &snpFile, SVParser &svFile);
  ~METHParser();

  void writeResult(PhasingResult phasingResult);
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
  // SNP map and iter
  std::map<int, SnpVariant> *currentVariants;
  std::map<int, SnpVariant>::iterator firstVariantIter;
  // SV map and iter
  std::map<int, std::map<int, bool>> *currentSV;
  std::vector<std::pair<int, int>>::iterator firstSVIter;
  std::map<std::string, std::vector<std::pair<int, int>>> SV_map;
  // mod map and iter
  std::map<int, std::map<std::string, ModVariant>> *currentMod;
  std::map<int, std::map<std::string, ModVariant>>::iterator firstModIter;
  void initSnpVariants(SnpParser &snpMap, const std::string &ref_string);
  void initSvVariants(SVParser &svFile);
  void initModVariants(METHParser &modFile);
  bool handleDeletionHomopolymerSNP(
      const bam1_t &aln, int ref_pos, int query_pos, int del_len,
      const std::string &ref_string,
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
               const std::string &ref_string, bool isONT, int svWindow,
               double svThreshold);
  void getClip(int pos, int clipFrontBack, int len, ClipCount &clipCount);
  bool isValidAlignment(const bam1_t &aln, int mappingQuality);
  void processBamFile(const std::string &bamFile, int lastSNPPos,
                      htsThreadPool &threadPool, const PhasingParameters &params,
                      std::vector<ReadVariant> &readVariantVec,
                      ClipCount &clipCount, const std::string &ref_string);

public:
  BamParser(std::string chrName, std::vector<std::string> inputBamFileVec,
            SnpParser &snpMap, SVParser &svFile, METHParser &modFile,
            const std::string &ref_string);
  ~BamParser();

  void direct_detect_alleles(int lastSNPPos, htsThreadPool &threadPool,
                             PhasingParameters params,
                             std::vector<ReadVariant> &readVariantVec,
                             ClipCount &clipCount,
                             const std::string &ref_string);
};

#endif
