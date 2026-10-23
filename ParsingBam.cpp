#include "ParsingBam.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>   // std::fixed / std::setprecision for INFO/PE|H1|H2
#include <sstream>
#include <string.h>

// vcf parser modify from
// http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html
// bam parser modify from
// https://github.com/whatshap/whatshap/blob/9882248c722b1020321fea6e9491c1cc5b75354b/whatshap/variants.py

// FASTA
FastaParser::FastaParser(std::string fastaFile,
                         std::vector<std::string> chrName,
                         std::vector<int> last_pos, int numThreads)
    : fastaFile(fastaFile), chrName(chrName), last_pos(last_pos) {
  // init map
  for (std::vector<std::string>::iterator iter = chrName.begin();
       iter != chrName.end(); iter++)
    chrString.insert(std::make_pair((*iter), ""));

  // load reference index
  faidx_t *fai = NULL;
  fai = fai_load(fastaFile.c_str());

  // iterating all chr
  // #pragma omp parallel for schedule(dynamic) num_threads(numThreads)
  for (std::vector<std::string>::iterator iter = chrName.begin();
       iter != chrName.end(); iter++) {

    int index = iter - chrName.begin();

    // Do not extract references without SNP coverage.
    if (last_pos.at(index) == -1) {
      chrString[(*iter)] = "";
      continue;
    }

    // ref_len is a return value that is length of retrun string
    int ref_len = 0;

    // read file; faidx_fetch_seq returns a malloc'd buffer owned by us
    char *seq = faidx_fetch_seq(fai, (*iter).c_str(), 0,
                                last_pos.at(index) + 5, &ref_len);
    if (ref_len == 0) {
      std::cout << "nothing in reference file \n";
    }

    // update map
    chrString[(*iter)] = seq ? std::string(seq, ref_len) : std::string();
    free(seq);
  }
  fai_destroy(fai);
}

FastaParser::~FastaParser() {}

void BaseVairantParser::readGzLines(
    const std::string &variantFile,
    std::function<void(const std::string &)> callback) {
  gzFile file = gzopen(variantFile.c_str(), "rb");
  if (!file) {
    std::cout << "Fail to open vcf: " << variantFile << "\n";
    return;
  }
  int buffer_size = 1048576; // 1M
  char *buffer = (char *)malloc(buffer_size);
  if (!buffer) {
    std::cerr << "Failed to allocate buffer\n";
    exit(EXIT_FAILURE);
  }
  char *offset = buffer;

  while (true) {
    int len = buffer_size - (offset - buffer);
    if (len == 0) {
      buffer_size *= 2; // Double the buffer size
      char *new_buffer = (char *)realloc(buffer, buffer_size);
      if (!new_buffer) {
        std::cerr << "Failed to allocate buffer\n";
        free(buffer);
        exit(EXIT_FAILURE);
      }
      buffer = new_buffer;
      offset =
          buffer +
          buffer_size /
              2; // Update the offset pointer to the end of the old buffer
      len = buffer_size - (offset - buffer);
    }

    len = gzread(file, offset, len);
    if (len == 0)
      break;
    if (len < 0) {
      int err;
      fprintf(stderr, "Error: %s.\n", gzerror(file, &err));
      exit(EXIT_FAILURE);
    }

    char *cur = buffer;
    char *end = offset + len;
    for (char *eol; (cur < end) && (eol = std::find(cur, end, '\n')) < end;
         cur = eol + 1) {
      callback(std::string(cur, eol));
    }
    // any trailing data in [eol, end) now is a partial line
    offset = std::copy(cur, end, buffer);
  }
  gzclose(file);
  free(buffer);
}

void BaseVairantParser::compressParser(std::string &variantFile) {
  if (variantFile == "") return;
  readGzLines(variantFile, [this](const std::string &line) {
    parserProcess(const_cast<std::string &>(line));
  });
}

void BaseVairantParser::unCompressParser(std::string &variantFile) {
  std::ifstream originVcf(variantFile);
  if (variantFile == "")
    return;
  if (!originVcf.is_open()) {
    std::cout << "Fail to open vcf: " << variantFile << "\n";
    exit(1);
  } else {
    std::string input;
    while (!originVcf.eof()) {
      std::getline(originVcf, input);
      parserProcess(input);
    }
  }
}

void BaseVairantParser::compressInput(std::string variantFile,
                                      std::string resultFile,
                                      PhasingResult &phasingResult) {
  if (variantFile == "") return;
  std::ofstream resultVcf(resultFile);
  if (!resultVcf.is_open()) {
    std::cout << "Fail to open write file: " << resultFile << "\n";
    return;
  }
  bool ps_def = false;
  readGzLines(variantFile, [&](const std::string &line) {
    writeLine(const_cast<std::string &>(line), ps_def, resultVcf, phasingResult);
  });
}

void BaseVairantParser::unCompressInput(std::string variantFile,
                                        std::string resultFile,
                                        PhasingResult &phasingResult) {
  std::ifstream originVcf(variantFile);
  std::ofstream resultVcf(resultFile);

  if (!resultVcf.is_open()) {
    std::cout << "Fail to open write file: " << resultFile << "\n";
  } else if (!originVcf.is_open()) {
    std::cout << "Fail to open vcf: " << variantFile << "\n";
  } else {
    bool ps_def = false;
    std::string input;
    while (!originVcf.eof()) {
      std::getline(originVcf, input);
      if (input != "") {
        writeLine(input, ps_def, resultVcf, phasingResult);
      }
    }
  }
}

void BaseVairantParser::dispatchWriteResult(const std::string &inputFile,
                                            const std::string &outputFile,
                                            PhasingResult &phasingResult) {
  if (inputFile.find("gz") != std::string::npos)
    compressInput(inputFile, outputFile, phasingResult);
  else if (inputFile.find("vcf") != std::string::npos)
    unCompressInput(inputFile, outputFile, phasingResult);
}

void BaseVairantParser::writeLine(std::string &input, bool &ps_def,
                                  std::ofstream &resultVcf,
                                  PhasingResult &phasingResult) {
  if (input.substr(0, 2) == "##") {
    // Drop PE/H1/H2 declarations carried over from an input that was
    // already phased by LongPhase; writePeInfoHeaders() emits the current
    // ones, so re-phasing never produces duplicate ##INFO lines.
    if (input.compare(0, 14, "##INFO=<ID=PE,") == 0 ||
        input.compare(0, 14, "##INFO=<ID=H1,") == 0 ||
        input.compare(0, 14, "##INFO=<ID=H2,") == 0)
      return;
    writeMetaHeader(input, ps_def, resultVcf);
  }
  else if (input.substr(0, 6) == "#CHROM" || input.substr(0, 6) == "#chrom")
    writeColumnHeader(input, ps_def, resultVcf);
  else
    writeDataLine(input, resultVcf, phasingResult);
}

void BaseVairantParser::writeMetaHeader(const std::string &input, bool &ps_def,
                                        std::ofstream &resultVcf) {
  if (input.substr(0, 16) == "##FORMAT=<ID=PS,")
    ps_def = true;
  resultVcf << input << "\n";
}

void BaseVairantParser::writeColumnHeader(const std::string &input,
                                          bool &ps_def,
                                          std::ofstream &resultVcf) {
  if (commandLine == false) {
    if (ps_def == false) {
      resultVcf << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description="
                   "\"Phase set identifier\">\n";
      ps_def = true;
    }
    // Emit PE/H1/H2 header lines here too as a safety net: if the input
    // VCF has no ##FILTER=<ID=PASS line (SnpParser hooks off that) or
    // no ##INFO block at all, they still get declared before #CHROM.
    writePeInfoHeaders(resultVcf);
    resultVcf << "##longphaseVersion=" << params->version << "\n";
    resultVcf << "##commandline=\"" << params->command << "\"\n";
    commandLine = true;
  }
  resultVcf << input << "\n";
}

void BaseVairantParser::writePeInfoHeaders(std::ofstream &resultVcf) {
  // Idempotent: writes the three ##INFO lines the first time we're
  // called, then never again for this parser instance. Every parser
  // (SNP, SV, MOD) shares this via inheritance.
  if (pe_def) return;
  resultVcf << "##INFO=<ID=PE,Number=1,Type=Float,Description=\"Phasing entropy in bits; 0 means either a unanimous vote or no incoming votes (block start), 1 means maximally ambiguous\">\n";
  resultVcf << "##INFO=<ID=H1,Number=1,Type=Float,Description=\"Weighted HP1 vote count\">\n";
  resultVcf << "##INFO=<ID=H2,Number=1,Type=Float,Description=\"Weighted HP2 vote count\">\n";
  pe_def = true;
}

namespace {
// Remove existing PE=, H1= and H2= entries from an INFO field, e.g. when
// re-phasing a VCF that LongPhase already annotated. Keys are matched
// exactly, so fields such as H1x= are kept. An emptied field becomes ".".
void removePeHapFromInfoField(std::string &infoField) {
  if (infoField == ".")
    return;
  std::string kept;
  bool removed = false;
  size_t start = 0;
  while (start <= infoField.size()) {
    size_t end = infoField.find(';', start);
    if (end == std::string::npos)
      end = infoField.size();
    std::string entry = infoField.substr(start, end - start);
    std::string key = entry.substr(0, entry.find('='));
    if (key == "PE" || key == "H1" || key == "H2") {
      removed = true;
    } else if (!entry.empty()) {
      if (!kept.empty())
        kept += ";";
      kept += entry;
    }
    start = end + 1;
  }
  // Leave fields without PE/H1/H2 untouched (unphased records stay as-is)
  if (removed)
    infoField = kept.empty() ? "." : kept;
}

// Format the PE/H1/H2 fragment and either replace a "." INFO field or
// append with a leading ';' to an existing INFO field. Any PE/H1/H2
// already present are removed first so each key appears only once.
// Shared by SNP, SV, and MOD writers so all three VCFs use the exact
// same formatting.
void appendPeHapInfoToInfoField(std::string &infoField,
                                const PhasingElement &pe) {
  removePeHapFromInfoField(infoField);
  std::ostringstream buf;
  buf << std::fixed << std::setprecision(3);
  buf << "PE=" << pe.entropy
      << ";H1=" << pe.h1
      << ";H2=" << pe.h2;
  // modcall ends INFO with ';', so drop trailing separators to avoid
  // writing an empty entry (";;PE=...")
  while (!infoField.empty() && infoField.back() == ';')
    infoField.pop_back();
  if (infoField.empty() || infoField == ".") {
    infoField = buf.str();
  } else {
    infoField += ";" + buf.str();
  }
}
} // namespace

// Move PS to the end of FORMAT in every sample column and set the PS of
// column sampleCol to value. The other samples keep the PS they had (or
// none), so phasing one sample leaves the others unchanged.
static void setSamplePs(std::vector<std::string> &fields, size_t sampleCol,
                        const std::string &value) {
  auto split = [](const std::string &s) {
    std::vector<std::string> v;
    size_t b = 0;
    while (true) {
      size_t e = s.find(':', b);
      v.push_back(s.substr(b, e == std::string::npos ? std::string::npos : e - b));
      if (e == std::string::npos)
        break;
      b = e + 1;
    }
    return v;
  };
  auto join = [](const std::vector<std::string> &v) {
    std::string s;
    for (size_t i = 0; i < v.size(); i++)
      s += (i ? ":" : "") + v[i];
    return s;
  };
  std::vector<std::string> keys = split(fields[8]);
  std::vector<std::string> oldPs(fields.size());
  for (size_t k = 0; k < keys.size(); k++) {
    if (keys[k] != "PS")
      continue;
    keys.erase(keys.begin() + k);
    for (size_t c = 9; c < fields.size(); c++) {
      std::vector<std::string> v = split(fields[c]);
      if (v.size() > k) {
        oldPs[c] = v[k];
        v.erase(v.begin() + k);
        fields[c] = join(v);
      }
    }
    break;
  }
  keys.push_back("PS");
  fields[8] = join(keys);
  for (size_t c = 9; c < fields.size(); c++) {
    const std::string &ps = (c == sampleCol) ? value : oldPs[c];
    if (ps.empty())
      continue; // no PS before: the trailing PS value stays omitted
    std::vector<std::string> v = split(fields[c]);
    while (v.size() + 1 < keys.size())
      v.push_back(".");
    v.push_back(ps);
    fields[c] = join(v);
  }
}

static bool isHetGt(const int *gt) {
  return (gt[0] == 2 && gt[1] == 4) || // 0/1
         (gt[0] == 4 && gt[1] == 2) || // 1/0
         (gt[0] == 2 && gt[1] == 5) || // 0|1
         (gt[0] == 4 && gt[1] == 3);   // 1|0
}

static void parseReadList(
    const std::string &info, const std::string &tag,
    const std::string &chr, int repPos, bool is_reverse, bool is_modify,
    std::map<std::string, std::map<int, std::map<std::string, ModVariant>>> &chrVariant) {
  size_t read_pos = info.find(tag);
  read_pos = info.find("=", read_pos) + 1;
  size_t next_field = info.find(";", read_pos);
  std::string totalRead = info.substr(read_pos, next_field - read_pos);
  std::stringstream ss(totalRead);
  std::string read;
  while (std::getline(ss, read, ',')) {
    ModVariant tmp;
    tmp.is_reverse = is_reverse;
    tmp.is_modify = is_modify;
    chrVariant[chr][repPos][read] = tmp;
  }
}

// SNP
SnpParser::SnpParser(PhasingParameters &in_params) {

  chrVariant = new std::map<std::string, std::map<int, SnpVariant>>;

  params = &in_params;

  // Initialize removed indels log file
  if (params->phaseIndel && params->indelQuality > 0) {
    std::string logFileName = params->resultPrefix + "_removed_indels.log";
    removedIndelsLog.open(logFileName.c_str());
    if (removedIndelsLog.is_open()) {
      removedIndelsLog << "#CHROM\tPOS\tREF\tALT\tQUAL\n";
    }
  }

  // open vcf file
  htsFile *inf = bcf_open(params->snpFile.c_str(), "r");
  // read header
  bcf_hdr_t *hdr = bcf_hdr_read(inf);
  // counters
  int nseq = 0;
  // report names of all the sequences in the VCF file
  const char **seqnames = NULL;
  // chromosome idx and name
  seqnames = bcf_hdr_seqnames(hdr, &nseq);
  // store chromosome
  for (int i = 0; i < nseq; i++) {
    // bcf_hdr_id2name is another way to get the name of a sequence
    chrName.push_back(seqnames[i]);
  }
  // the sample to phase (--sample), by default the first one
  if (bcf_hdr_nsamples(hdr) == 0) {
    std::cerr << "error: " << params->snpFile << " has no sample column\n";
    exit(1);
  }
  std::string sampleName = params->sample.empty() ? std::string(hdr->samples[0]) : params->sample;
  int sampleIdx = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, sampleName.c_str());
  if (sampleIdx < 0) {
    std::cerr << "error: sample " << sampleName << " is not in " << params->snpFile
              << ". samples:";
    for (int i = 0; i < bcf_hdr_nsamples(hdr); i++)
      std::cerr << " " << hdr->samples[i];
    std::cerr << "\n";
    exit(1);
  }
  // column of this sample in the VCF text, used when writing the result
  sampleCol = 9 + sampleIdx;
  // limit the VCF data to this sample
  int is_file = bcf_hdr_set_samples(hdr, sampleName.c_str(), 0);
  if (is_file != 0) {
    std::cout << "error or a positive integer if the list contains samples not "
                 "present in the VCF header\n";
  }

  // struct for storing each record
  bcf1_t *rec = bcf_init();
  int ngt_arr = 0;
  int ngt = 0;
  int *gt = NULL;

  // loop vcf line
  while (bcf_read(inf, hdr, rec) == 0) {
    // snp
    if (bcf_is_snp(rec)) {
      ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);

      if (ngt < 0) {
        std::cerr << "pos " << rec->pos << " missing GT value" << "\n";
        exit(1);
      }

      // just phase hetero SNP
      if (isHetGt(gt)) {

        // get chromosome string
        std::string chr = seqnames[rec->rid];
        // position is 0-base
        int variantPos = rec->pos;
        // get r alleles
        SnpVariant tmp;
        tmp.Ref = rec->d.allele[0];
        tmp.Alt = rec->d.allele[1];

        // prevent the MAVs calling error which makes the GT=0/1
        if (rec->d.allele[1][2] != '\0') {
          continue;
        }

        // record
        (*chrVariant)[chr][variantPos] = tmp;
      }
    }
    // indel
    else if (params->phaseIndel) {
      ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);

      if (ngt < 0) {
        std::cerr << "pos " << rec->pos << " missing GT value" << "\n";
        exit(1);
      }

      if (isHetGt(gt)) {

        // get chromosome string
        std::string chr = seqnames[rec->rid];
        // position is 0-base
        int variantPos = rec->pos;
        // get r alleles
        SnpVariant tmp;
        tmp.Ref = rec->d.allele[0];
        tmp.Alt = rec->d.allele[1];

        float qual = rec->qual;
        if (std::isnan(qual)) {
          qual = 0.0;
        }
        if (params->indelQuality > 0 && qual < params->indelQuality) {
          if (removedIndelsLog.is_open()) {
            removedIndelsLog << chr << "\t" << (variantPos + 1) << "\t"
                             << tmp.Ref << "\t" << tmp.Alt << "\t"
                             << (std::isnan(rec->qual) ? "." : std::to_string(rec->qual))
                             << "\n";
          }
          // Record the position of filtered indels (0-based)
          filteredIndelPositions[chr].insert(variantPos);
          continue;
        }

        // prevent the MAVs calling error which makes the GT=0/1
        if (rec->d.allele[1][tmp.Alt.size() + 1] != '\0') {
          continue;
        }

        // record
        (*chrVariant)[chr][variantPos] = tmp;
      }
    }
  }
}

SnpParser::~SnpParser() {
  if (removedIndelsLog.is_open()) {
    removedIndelsLog.close();
  }
  delete chrVariant;
}

std::map<int, SnpVariant> SnpParser::getVariants(std::string chrName) {
  std::map<int, SnpVariant> targetVariants;
  std::map<std::string, std::map<int, SnpVariant>>::iterator chrIter =
      chrVariant->find(chrName);

  if (chrIter != chrVariant->end())
    targetVariants = (*chrIter).second;

  return targetVariants;
}

std::map<int, SnpVariant>
SnpParser::getVariants_markindel(std::string chrName, const std::string &ref) {
  std::map<int, SnpVariant> targetVariants;
  std::map<std::string, std::map<int, SnpVariant>>::iterator chrIter =
      chrVariant->find(chrName);

  // Mark the indel which lies in the tandem repeat
  if (chrIter != chrVariant->end()) {
    for (auto innerIter = chrIter->second.begin();
         innerIter != chrIter->second.end(); innerIter++) {
      innerIter->second.is_danger = isDangerIndel(
          innerIter->first, ref, innerIter->second.Ref.length(),
          innerIter->second.Alt.length());
    }

    targetVariants = (*chrIter).second;
  }

  return targetVariants;
}

std::vector<std::string> SnpParser::getChrVec() { return chrName; }


int SnpParser::getLastSNP(std::string chrName) {
  std::map<std::string, std::map<int, SnpVariant>>::iterator chrVariantIter =
      chrVariant->find(chrName);
  // this chromosome not exist in this file.
  if (chrVariantIter == chrVariant->end())
    return -1;
  // get last SNP
  std::map<int, SnpVariant>::reverse_iterator lastVariantIter =
      (*chrVariantIter).second.rbegin();
  // there are no SNPs on this chromosome
  if (lastVariantIter == (*chrVariantIter).second.rend())
    return -1;
  return (*lastVariantIter).first;
}

void SnpParser::writeResult(PhasingResult &phasingResult,
                            const std::string &outPrefix) {
  dispatchWriteResult(params->snpFile, outPrefix + ".vcf",
                      phasingResult);
}

void SnpParser::parserProcess(std::string &input) {}

void SnpParser::writeMetaHeader(const std::string &input, bool &ps_def,
                               std::ofstream &resultVcf) {
  if (input.substr(0, 16) == "##FORMAT=<ID=PS,") {
    ps_def = true;
  }
  if (input.substr(0, 17) == "##FILTER=<ID=PASS") {
    resultVcf << input << "\n";
    if (params->phaseIndel && params->indelQuality > 0) {
      resultVcf << "##FILTER=<ID=INDEL_QUAL_FILTERED,Description=\"Indel "
                   "filtered due to QUAL below threshold ("
                << params->indelQuality << ")\">\n";
    }
    // Emit PE/H1/H2 INFO headers right after ##FILTER=<ID=PASS so they
    // sit in a natural place near the other INFO declarations. Safe if
    // called more than once (writePeInfoHeaders is idempotent).
    writePeInfoHeaders(resultVcf);
  } else {
    resultVcf << input << "\n";
  }
}


void SnpParser::writeDataLine(const std::string &input,
                              std::ofstream &resultVcf,
                              PhasingResult &phasingResult) {
  std::istringstream iss(input);
  std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),
                                  std::istream_iterator<std::string>());

  if (fields.size() == 0)
    return;

  // the phased sample's column (--sample); other sample columns are kept
  const size_t sc = sampleCol;
  if (fields.size() <= sc) {
    resultVcf << input << "\n";
    return;
  }

  int pos = std::stoi(fields[1]);
  int posIdx = pos - 1;
  std::string key = fields[0] + "_" + std::to_string(posIdx);

  PhasingResult::iterator psElementIter = phasingResult.find(key);

  // reset GT flag
  if (fields[8].find("GT") != std::string::npos) {
    int modify_start = findFieldValueStart(fields[8], fields[sc], "GT");
    if (fields[sc][modify_start + 1] == '|') {
      // direct modify GT value
      if (fields[sc][modify_start] > fields[sc][modify_start + 2]) {
        fields[sc][modify_start + 1] = fields[sc][modify_start];
        fields[sc][modify_start] = fields[sc][modify_start + 2];
        fields[sc][modify_start + 2] = fields[sc][modify_start + 1];
      }
      fields[sc][modify_start + 1] = '/';
    }
  }

  // Check if the variant is extracted from this VCF
  auto posIter = (*chrVariant)[fields[0]].find(posIdx);

  // Check if this indel was filtered out due to quality
  bool isFilteredIndel = false;
  if (params->phaseIndel && params->indelQuality > 0) {
    auto filteredIter = filteredIndelPositions.find(fields[0]);
    if (filteredIter != filteredIndelPositions.end()) {
      if (filteredIter->second.find(posIdx) != filteredIter->second.end()) {
        isFilteredIndel = true;
      }
    }
  }

  // this pos is phase; only the record that was phased (same REF and ALT)
  // gets the phase, not other records at the same position (issue #133)
  if (psElementIter != phasingResult.end() &&
      posIter != (*chrVariant)[fields[0]].end() &&
      posIter->second.Ref == fields[3] && posIter->second.Alt == fields[4]) {
    // add PS flag and value
    setSamplePs(fields, sc, std::to_string((*psElementIter).second.block));

    int modify_start = findFieldValueStart(fields[8], fields[sc], "GT");
    // direct modify GT value
    fields[sc][modify_start] = (*psElementIter).second.RAstatus[0];
    fields[sc][modify_start + 1] = '|';
    fields[sc][modify_start + 2] = (*psElementIter).second.RAstatus[2];

    // Attach PE / H1 / H2 to INFO for phased positions only.
    // Rule: if INFO was ".", replace it; otherwise append with ';'.
    appendPeHapInfoToInfoField(fields[7], (*psElementIter).second);
  }
  // this pos has not been phased
  else {
    // drop stale PE/H1/H2 left over from a previously phased input
    removePeHapFromInfoField(fields[7]);
    // add PS flag and value
    setSamplePs(fields, sc, ".");
  }

  // Add FILTER tag for filtered indels
  if (isFilteredIndel) {
    // Overwrite to INDEL_QUAL_FILTERED, do not keep the original tag
    fields[6] = "INDEL_QUAL_FILTERED";
  }

  for (std::vector<std::string>::iterator fieldIter = fields.begin();
       fieldIter != fields.end(); ++fieldIter) {
    if (fieldIter != fields.begin())
      resultVcf << "\t";
    resultVcf << (*fieldIter);
  }
  resultVcf << "\n";
}

bool SnpParser::findSNP(std::string chr, int position) {
  std::map<std::string, std::map<int, SnpVariant>>::iterator chrIter =
      chrVariant->find(chr);
  // empty chromosome
  if (chrIter == chrVariant->end())
    return false;

  std::map<int, SnpVariant>::iterator posIter =
      (*chrVariant)[chr].find(position);
  // empty position
  if (posIter == (*chrVariant)[chr].end())
    return false;

  return true;
}

void SnpParser::filterSNP(std::string chr,
                          std::vector<ReadVariant> &readVariantVec,
                          std::string &chr_reference) {

  // pos, <allele, <strand, True>
  std::map<int, std::map<int, std::map<int, bool>>> posAlleleStrand;
  std::map<int, bool> methylation;

  /*
  // iter all variant, record the strand contained in each SNP
  for( auto readSNPVecIter : readVariantVec ){
      // tag allele on forward or reverse strand
      for(auto variantIter : readSNPVecIter.variantVec ){
          posAlleleStrand[variantIter.position][variantIter.allele][readSNPVecIter.is_reverse]
  = true;
      }
  }

  // iter all SNP, both alleles that require SNP need to appear in the two
  strand for(auto pos: posAlleleStrand){
      // this position contain two allele, REF allele appear in the two strand,
  ALT allele appear in the two strand if(pos.second.size() == 2 &&
  pos.second[0].size() == 2 && pos.second[1].size() == 2 ){
          // high confident SNP
      }
      else{
          //methylation[pos.first] = true;
      }
  }
  */

  // Filter SNPs that are not easy to phasing due to homopolymer
  // get variant list
  std::map<std::string, std::map<int, SnpVariant>>::iterator chrIter =
      chrVariant->find(chr);
  std::map<int, bool> errorProneSNP;

  if (chrIter != chrVariant->end()) {
    std::map<int, int> consecutiveAllele;
    // iter all SNP and tag homopolymer
    for (std::map<int, SnpVariant>::iterator posIter =
             (*chrIter).second.begin();
         posIter != (*chrIter).second.end(); posIter++) {
      consecutiveAllele[(*posIter).first] =
          homopolymerLength((*posIter).first, chr_reference);
    }
    std::map<int, SnpVariant>::iterator currSNPIter = (*chrIter).second.begin();
    std::map<int, SnpVariant>::iterator nextSNPIter = std::next(currSNPIter, 1);
    // check whether each SNP pair falls in an area that is not easy to phasing
    while (currSNPIter != (*chrIter).second.end() &&
           nextSNPIter != (*chrIter).second.end()) {
      int currPos = (*currSNPIter).first;
      int nextPos = (*nextSNPIter).first;
      // filter one of SNP if this SNP pair falls in homopolymer and distance<=2
      if (consecutiveAllele[currPos] >= 3 && consecutiveAllele[nextPos] >= 3 &&
          std::abs(currPos - nextPos) <= 2) {
        errorProneSNP[nextPos] = true;
        nextSNPIter = (*chrIter).second.erase(nextSNPIter);
        continue;
      }

      currSNPIter++;
      nextSNPIter++;
    }
  }

  // iter all reads
  for (std::vector<ReadVariant>::iterator readSNPVecIter =
           readVariantVec.begin();
       readSNPVecIter != readVariantVec.end(); readSNPVecIter++) {
    // iter all SNPs in this read
    for (std::vector<Variant>::iterator variantIter =
             (*readSNPVecIter).variantVec.begin();
         variantIter != (*readSNPVecIter).variantVec.end();) {
      std::map<int, bool>::iterator delSNPIter =
          methylation.find((*variantIter).position);
      std::map<int, bool>::iterator homoIter =
          errorProneSNP.find((*variantIter).position);

      if (delSNPIter != methylation.end()) {
        variantIter = (*readSNPVecIter).variantVec.erase(variantIter);
      } else if (homoIter != errorProneSNP.end()) {
        variantIter = (*readSNPVecIter).variantVec.erase(variantIter);
      } else {
        variantIter++;
      }
    }
  }
}

// value of an INFO key, matching the whole key between ';'
static bool getInfoValue(const std::string &info, const std::string &key,
                         std::string &value) {
  size_t start = 0;
  while (start < info.size()) {
    size_t end = info.find(';', start);
    if (end == std::string::npos)
      end = info.size();
    if (end - start > key.size() && info.compare(start, key.size(), key) == 0 &&
        info[start + key.size()] == '=') {
      value = info.substr(start + key.size() + 1, end - start - key.size() - 1);
      return true;
    }
    start = end + 1;
  }
  return false;
}

// SV
SVParser::SVParser(PhasingParameters &in_params, SnpParser &in_snpFile) {
  params = &in_params;
  snpFile = &in_snpFile;

  chrVariant = new std::map<std::string, std::map<int, std::map<int, bool>>>;

  if (params->svFile.find("gz") != std::string::npos) {
    // .vcf.gz
    compressParser(params->svFile);
  } else if (params->svFile.find("vcf") != std::string::npos) {
    // .vcf
    unCompressParser(params->svFile);
  }

  for (auto &chrRecords : chrDeletionRecords)
    std::sort(chrRecords.second.begin(), chrRecords.second.end());

  // erase SV pos if this pos appear two or more times
  for (std::map<std::string, std::map<int, bool>>::iterator chrIter =
           posDuplicate.begin();
       chrIter != posDuplicate.end(); chrIter++) {
    for (std::map<int, bool>::iterator posIter = (*chrIter).second.begin();
         posIter != (*chrIter).second.end(); posIter++) {
      if ((*posIter).second == true) {
        std::map<int, std::map<int, bool>>::iterator erasePosIter =
            (*chrVariant)[(*chrIter).first].find((*posIter).first);
        if (erasePosIter != (*chrVariant)[(*chrIter).first].end()) {
          (*chrVariant)[(*chrIter).first].erase(erasePosIter);
        }
      }
    }
  }

}

SVParser::~SVParser() { delete chrVariant; }

void SVParser::parserProcess(std::string &input) {
  if (input.substr(0, 1) == "#")
    return;

  std::istringstream iss(input);
  std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),
                                  std::istream_iterator<std::string>());

  if (fields.size() == 0)
    return;
  int pos = std::stoi(fields[1]) - 1;
  int start = std::stoi(fields[1]);
  std::string chr = fields[0];

  int modify_start = findFieldValueStart(fields[8], fields[9], "GT");
  bool filter = false;

  // homo GT
  if (fields[9][modify_start] == fields[9][modify_start + 2]) {
    filter = true;
  }
  // conflict pos with SNP
  if ((*snpFile).findSNP(chr, pos)) {
    filter = true;
  }

  std::map<int, bool>::iterator posIter = posDuplicate[chr].find(pos);
  // conflict pos with SV
  if (posIter == posDuplicate[chr].end())
    posDuplicate[chr][pos] = false;
  else {
    posDuplicate[chr][pos] = true;
    filter = true;
  }

  // every deletion record, whatever its genotype and whether phase uses it,
  // so that a read's split junction is given to the nearest of them
  // (issue #126)
  std::string recordSvlen;
  if (getInfoValue(fields[7], "SVLEN", recordSvlen)) {
    char *svlenEnd;
    long long recordLength = std::strtoll(recordSvlen.c_str(), &svlenEnd, 10);
    std::string recordType;
    bool recordIsDeletion = getInfoValue(fields[7], "SVTYPE", recordType)
                                ? recordType == "DEL"
                                : recordLength < 0;
    if (svlenEnd != recordSvlen.c_str() && recordIsDeletion &&
        recordLength != 0 && recordLength > -INT_MAX && recordLength < INT_MAX)
      chrDeletionRecords[chr].emplace_back(start, (int)std::llabs(recordLength));
  }

  if (filter) {
    return;
  }

  // get read INFO
  std::string info = fields[7];
  size_t svlenPos = info.find("SVLEN=");
  if (svlenPos != std::string::npos) {
    svlenPos += 6;
    size_t semiPos = info.find(';', svlenPos);
    int svlen = std::stoi(info.substr(svlenPos, semiPos - svlenPos));
  (*chrVariant)[chr][start][svlen] = true;

    // SV type and the reads the SV caller found supporting the SV (issue #126)
    SVReadSupport &support = chrSVSupport[chr][start];
    std::string svtype;
    if (getInfoValue(info, "SVTYPE", svtype))
      support.isDeletion = (svtype == "DEL");
    else
      support.isDeletion = (svlen < 0);
    std::string rnames;
    if (getInfoValue(info, "RNAMES", rnames)) {
      std::istringstream rnamesStream(rnames);
      std::string readName;
      while (std::getline(rnamesStream, readName, ','))
        if (!readName.empty() && readName != ".")
          support.readNames.insert(readName);
      support.hasReadNames = !support.readNames.empty();
    }
  }
}

std::map<int, std::map<int, bool>> SVParser::getVariants(std::string chrName) {
  std::map<int, std::map<int, bool>> targetVariants;
  std::map<std::string, std::map<int, std::map<int, bool>>>::iterator chrIter =
      chrVariant->find(chrName);

  if (chrIter != chrVariant->end())
    targetVariants = (*chrIter).second;

  return targetVariants;
}

const std::map<int, SVReadSupport> *
SVParser::getReadSupport(const std::string &chrName) const {
  std::map<std::string, std::map<int, SVReadSupport>>::const_iterator chrIter =
      chrSVSupport.find(chrName);
  if (chrIter == chrSVSupport.end())
    return nullptr;
  return &(*chrIter).second;
}

const std::vector<std::pair<int, int>> *
SVParser::getDeletionRecords(const std::string &chrName) const {
  std::map<std::string, std::vector<std::pair<int, int>>>::const_iterator
      chrIter = chrDeletionRecords.find(chrName);
  if (chrIter == chrDeletionRecords.end())
    return nullptr;
  return &(*chrIter).second;
}

void SVParser::writeResult(PhasingResult &phasingResult,
                           const std::string &outPrefix) {
  dispatchWriteResult(params->svFile, outPrefix + "_SV.vcf",
                      phasingResult);
}

void SVParser::writeMetaHeader(const std::string &input, bool &ps_def,
                               std::ofstream &resultVcf) {
  if (input.substr(0, 16) == "##FORMAT=<ID=PS,") {
    ps_def = true;
  }
  // SV VCFs vary widely (Sniffles / cuteSV / DELLY etc.) and may or may
  // not have ##FILTER=<ID=PASS. To be robust, hook off any ##INFO line
  // and slip the PE/H1/H2 declarations in alongside them. If no ##INFO
  // lines exist at all, writeColumnHeader's safety net will still emit
  // them before #CHROM.
  resultVcf << input << "\n";
  if (input.substr(0, 7) == "##INFO=") {
    writePeInfoHeaders(resultVcf);
  }
}

void SVParser::writeDataLine(const std::string &input,
                               std::ofstream &resultVcf,
                               PhasingResult &phasingResult) {
  std::istringstream iss(input);
  std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),
                                  std::istream_iterator<std::string>());

  if (fields.size() == 0)
    return;

  int pos = std::stoi(fields[1]);
  int posIdx = pos - 1;
  std::string key = fields[0] + "_" + std::to_string(posIdx);

  PhasingResult::iterator psElementIter = phasingResult.find(key);

  // PS flag already exist, erase PS info
  if (fields[8].find("PS") != std::string::npos) {
    std::string format_before_ps_erase = fields[8];
    int ps_pos = fields[8].find("PS");

    // erase PS flag
    if (fields[8].find(":", ps_pos + 1) != std::string::npos) {
      fields[8].erase(ps_pos, 3);
    } else {
      fields[8].erase(ps_pos - 1, 3);
    }

    int ps_start = findFieldValueStart(format_before_ps_erase, fields[9], "PS");

    // erase PS value
    if (fields[9].find(":", ps_start + 1) != std::string::npos) {
      int ps_end_pos = fields[9].find(":", ps_start + 1);
      fields[9].erase(ps_start, ps_end_pos - ps_start + 1);
    } else {
      fields[9].erase(ps_start - 1, fields[9].length() - ps_start + 1);
    }
  }
  // reset GT flag
  if (fields[8].find("GT") != std::string::npos) {
    int modify_start = findFieldValueStart(fields[8], fields[9], "GT");
    if (fields[9][modify_start + 1] == '|') {
      // direct modify GT value
      if (fields[9][modify_start] > fields[9][modify_start + 2]) {
        fields[9][modify_start + 1] = fields[9][modify_start];
        fields[9][modify_start] = fields[9][modify_start + 2];
        fields[9][modify_start + 2] = fields[9][modify_start + 1];
      }
      fields[9][modify_start + 1] = '/';
    }
  }

  // Check if the variant is extracted from this VCF
  auto posIter = (*chrVariant)[fields[0]].find(posIdx + 1);

  // this pos is phase and exist in map
  if (psElementIter != phasingResult.end() &&
      posIter != (*chrVariant)[fields[0]].end()) {

    // add PS flag and value
    fields[8] = fields[8] + ":PS";
    fields[9] =
        fields[9] + ":" + std::to_string((*psElementIter).second.block);

    int modify_start = findFieldValueStart(fields[8], fields[9], "GT");
    // direct modify GT value
    fields[9][modify_start] = (*psElementIter).second.RAstatus[0];
    fields[9][modify_start + 1] = '|';
    fields[9][modify_start + 2] = (*psElementIter).second.RAstatus[2];

    // Attach PE / H1 / H2 to INFO for phased positions only.
    // Rule: if INFO was ".", replace it; otherwise append with ';'.
    appendPeHapInfoToInfoField(fields[7], (*psElementIter).second);
  }
  // this pos has not been phased
  else {
    // drop stale PE/H1/H2 left over from a previously phased input
    removePeHapFromInfoField(fields[7]);
    // add PS flag and value
    fields[8] = fields[8] + ":PS";
    fields[9] = fields[9] + ":.";
  }
  for (std::vector<std::string>::iterator fieldIter = fields.begin();
       fieldIter != fields.end(); ++fieldIter) {
    if (fieldIter != fields.begin())
      resultVcf << "\t";
    resultVcf << (*fieldIter);
  }
  resultVcf << "\n";
}

bool SVParser::findSV(std::string chr, int position) {
  std::map<std::string, std::map<int, std::map<int, bool>>>::iterator chrIter =
      chrVariant->find(chr);
  // empty chromosome
  if (chrIter == chrVariant->end())
    return false;

  std::map<int, std::map<int, bool>>::iterator posIter =
      (*chrVariant)[chr].find(position);
  // empty position
  if (posIter == (*chrVariant)[chr].end())
    return false;

  return true;
}
BamParser::BamParser(std::string inputChrName,
                     std::vector<std::string> inputBamFileVec,
                     SnpParser &snpMap, SVParser &svFile, METHParser &modFile,
                     const std::string &ref_string)
    : chrName(inputChrName), BamFileVec(inputBamFileVec) {

  currentVariants = new std::map<int, SnpVariant>;
  currentSV = new std::map<int, std::map<int, bool>>;
  currentMod = new std::map<int, std::map<std::string, ModVariant>>;

  initSnpVariants(snpMap, ref_string);
  initSvVariants(svFile);
  initModVariants(modFile);
}

BamParser::~BamParser() {
  delete currentVariants;
  delete currentSV;
  delete currentMod;
}

void BamParser::initSnpVariants(SnpParser &snpMap, const std::string &ref_string) {
  // use chromosome to find recorded snp map
  (*currentVariants) = snpMap.getVariants_markindel(chrName, ref_string);
  // set skip variant start iterator
  firstVariantIter = currentVariants->begin();
  if (firstVariantIter == currentVariants->end()) {
    std::cerr << "error chromosome name or empty map: " << chrName << "\n";
    exit(1);
  }
}

void BamParser::initSvVariants(SVParser &svFile) {
  // set current chromosome SV map
  (*currentSV) = svFile.getVariants(chrName);
  for (const auto &chrPair : *currentSV) {
    int start = chrPair.first;
    for (const auto &endPair : chrPair.second) {
      int svlen = endPair.first;
      SV_map[chrName].emplace_back(start, svlen);
    }
  }
  firstSVIter = SV_map[chrName].begin();
  currentSVSupport = svFile.getReadSupport(chrName);
  currentDeletionRecords = svFile.getDeletionRecords(chrName);
}

void BamParser::initModVariants(METHParser &modFile) {
  // set current chromosome MOD map
  (*currentMod) = modFile.getVariants(chrName);
  firstModIter = currentMod->begin();
}

bool BamParser::isValidAlignment(const bam1_t &aln, int mappingQuality) {
  // Intentionally different from modcall/haplotag filters: parameterized qual
  // threshold; supplementary (0x800) alignments are kept for phasing evidence.
  int flag = aln.core.flag;
  return !(aln.core.qual < mappingQuality // mapping quality
           || (flag & 0x4) != 0           // read unmapped
           || (flag & 0x100) != 0         // secondary alignment. repeat.
                                          // A secondary alignment occurs when a
                                          // given read could align reasonably
                                          // well to more than one place.
           || (flag & 0x400) != 0         // duplicate
           // (flag & 0x800) != 0         // supplementary alignment
           // A chimeric alignment is represented as a set of linear alignments
           // that do not have large overlaps.
  );
}

void BamParser::direct_detect_alleles(int lastSNPPos, htsThreadPool &threadPool,
                                      PhasingParameters params,
                                      std::vector<ReadVariant> &readVariantVec,
                                      ClipCount &clipCount,
                                      const std::string &ref_string) {

  // record SNP start iter
  std::map<int, SnpVariant>::iterator tmpFirstVariantIter = firstVariantIter;
  // record SV start iter
  std::vector<std::pair<int, int>>::iterator tmpFirstSVIter = firstSVIter;
  // record MOD start iter
  std::map<int, std::map<std::string, ModVariant>>::iterator tmpFirstModIter =
      firstModIter;

  for (auto bamFile : BamFileVec) {

    firstVariantIter = tmpFirstVariantIter;
    firstSVIter = tmpFirstSVIter;
    firstModIter = tmpFirstModIter;

    processBamFile(bamFile, lastSNPPos, threadPool, params, readVariantVec,
                   clipCount, ref_string);
  }
}

void BamParser::processBamFile(const std::string &bamFile, int lastSNPPos,
                                htsThreadPool &threadPool,
                                const PhasingParameters &params,
                                std::vector<ReadVariant> &readVariantVec,
                                ClipCount &clipCount,
                                const std::string &ref_string) {
  // open bam file
  samFile *fp_in = hts_open(bamFile.c_str(), "r");
  // load reference file
  hts_set_fai_filename(fp_in, params.fastaFile.c_str());
  // read header
  bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
  // initialize an alignment
  bam1_t *aln = bam_init1();
  hts_idx_t *idx = NULL;

  if ((idx = sam_index_load(fp_in, bamFile.c_str())) == 0) {
    std::cout << "ERROR: Cannot open index for bam file\n";
    exit(1);
  }

  std::string range = chrName + ":1-" + std::to_string(lastSNPPos);
  hts_itr_t *iter = sam_itr_querys(idx, bamHdr, range.c_str());

  hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &threadPool);
  int result;
  while ((result = sam_itr_multi_next(fp_in, iter, aln)) >= 0) {
    if (!isValidAlignment(*aln, params.mappingQuality)) {
      continue;
    }

      get_snp(*bamHdr, *aln, readVariantVec, clipCount, ref_string,
            params.isONT, params.svWindow, params.svThreshold,
            params.mappingQuality);
  }
  hts_idx_destroy(idx);
  bam_hdr_destroy(bamHdr);
  bam_destroy1(aln);
  sam_close(fp_in);
}

// Returns true if get_snp() should return immediately (read bounds exceeded).
bool BamParser::handleDeletionHomopolymerSNP(
    const bam1_t &aln, int ref_pos, int query_pos, int del_len,
    const std::string &ref_string,
    std::map<int, SnpVariant>::iterator &currentVariantIter,
    ReadVariant &tmpReadResult) {

  if (ref_pos + del_len + 1 == (*currentVariantIter).first) {
    return false;
  } else if ((*currentVariantIter).first >= ref_pos &&
             (*currentVariantIter).first < ref_pos + del_len) {
    // check snp in homopolymer
    if (homopolymerLength((*currentVariantIter).first, ref_string) >= 3) {
      int refAlleleLen = (*currentVariantIter).second.Ref.length();
      int altAlleleLen = (*currentVariantIter).second.Alt.length();
      int base_q = 0;

      if (query_pos + 1 > aln.core.l_qseq) {
        return true;
      }

      int allele = -1;
      // SNP
      if (refAlleleLen == 1 && altAlleleLen == 1) {
        // get the next match
        char base = seq_nt16_str[bam_seqi(bam_get_seq(&aln), query_pos)];
        if (base == (*currentVariantIter).second.Ref[0]) {
          allele = 0;
        } else if (base == (*currentVariantIter).second.Alt[0]) {
          allele = 1;
        }
        base_q = bam_get_qual(&aln)[query_pos];
      }
      // the read deletion contain VCF's deletion
      else if (refAlleleLen != 1 && altAlleleLen == 1) {
        allele = 1;
        // using this quality to identify indel
        base_q = -4;
      }

      if (allele != -1) {
        tmpReadResult.variantVec.push_back(
            Variant((*currentVariantIter).first, allele, base_q,
                    VariantType::SNP));
        currentVariantIter++;
      }
    }
  }
  return false;
}

void BamParser::handleSNPVariant(
    const bam1_t &aln, int cigar_idx, int cigar_oplen,
    int ref_pos, int query_pos,
    std::map<int, SnpVariant>::iterator &currentVariantIter,
    int &variantPos, ReadVariant &tmpReadResult) {

  const uint32_t *cigar = bam_get_cigar(&aln);
  int aln_core_n_cigar = aln.core.n_cigar;
  int refAlleleLen = (*currentVariantIter).second.Ref.length();
  int altAlleleLen = (*currentVariantIter).second.Alt.length();
  int offset = variantPos - ref_pos;
  int base_q = 0;
  int allele = -1;

  // SNP
  if (refAlleleLen == 1 && altAlleleLen == 1) {
    char base =
        seq_nt16_str[bam_seqi(bam_get_seq(&aln), query_pos + offset)];
    if (base == (*currentVariantIter).second.Ref[0])
      allele = 0;
    else if (base == (*currentVariantIter).second.Alt[0])
      allele = 1;
    base_q = bam_get_qual(&aln)[query_pos + offset];
  }
  // insertion
  else if (refAlleleLen == 1 && altAlleleLen != 1 &&
           cigar_idx + 1 < aln_core_n_cigar) {
    if (ref_pos + cigar_oplen - 1 == variantPos &&
        bam_cigar_op(cigar[cigar_idx + 1]) == 1) {
      allele = 1;
    } else {
      allele = 0;
    }
    // using this quality to identify indel
    base_q = -4;
    // using this quality to identify danger indel
    if ((*currentVariantIter).second.is_danger) {
      base_q = -5;
    }
  }
  // deletion
  else if (refAlleleLen != 1 && altAlleleLen == 1 &&
           cigar_idx + 1 < aln_core_n_cigar) {
    if (ref_pos + cigar_oplen - 1 == variantPos &&
        bam_cigar_op(cigar[cigar_idx + 1]) == 2) {
      allele = 1;
    } else {
      allele = 0;
    }
    // using this quality to identify indel
    base_q = -4;
    // using this quality to identify danger indel
    if ((*currentVariantIter).second.is_danger) {
      base_q = -5;
    }
  }

  if (allele != -1) {
    // record snp result
    tmpReadResult.variantVec.push_back(
        Variant(variantPos, allele, base_q, VariantType::SNP));
  }
  currentVariantIter++;
  if (currentVariantIter != currentVariants->end()) {
    variantPos = (*currentVariantIter).first;
  } else {
    variantPos = INT_MAX;
  }
}

void BamParser::handleModVariant(
    const bam1_t &aln, int variantPos,
    std::map<int, std::map<std::string, ModVariant>>::iterator &currentModIter,
    int &modPos, ReadVariant &tmpReadResult) {

  std::map<std::string, ModVariant>::iterator readIter =
      (*currentModIter).second.find(bam_get_qname(&aln));
  if (readIter != (*currentModIter).second.end() &&
      modPos <= variantPos) {
    // check variant strand in vcf file is same as bam file
    if ((*readIter).second.is_reverse == bam_is_rev(&aln)) {
      // -2 : forward strand
      // -3 : reverse strand
      int strand = (bam_is_rev(&aln) ? -3 : -2);
      int allele = ((*readIter).second.is_modify ? 0 : 1);
      // push mod into result vector
      tmpReadResult.variantVec.push_back(
          Variant(modPos, allele, strand, VariantType::MOD));
    }
  }
  currentModIter++;
  if (currentModIter != currentMod->end()) {
    modPos = (*currentModIter).first;
  } else {
    modPos = INT_MAX;
  }
}

// SV evidence from split alignments (issue #126): minimap2 reports a large
// deletion as a primary and a supplementary alignment (SA tag) rather than as a
// CIGAR D, so the alignment that reaches the SV ends, clipped, at its breakpoint.
// Breakpoint tolerance (bp) of the split-alignment junction; the aligner may
// extend an alignment over microhomology, which is added on top of it.
static const int SV_JUNCTION_TOLERANCE = 50;
// largest query overlap of two segments of a split read accepted at a junction
static const int SV_MAX_QUERY_OVERLAP = 1000;
// An alignment shows the reference at an SV if it is aligned at least this
// many bp on both sides of the SV position; such alignments are judged by the
// CIGAR alone, as before.
static const int SV_REF_FLANK = 100;
// a clip longer than this marks a break in the alignment (as in getClip)
static const int SV_MIN_BREAK_CLIP = 5;
// a split junction whose breakpoints both lie within this distance of the
// SV's may be the same deletion called at a slightly different place, and an
// alignment clipped within this distance (or |SVLEN|) past the SV may be the
// upstream part of a carrier whose alignment ran on over a repeat
static const int SV_MAX_JUNCTION_SHIFT = 1000;
// For a record without RNAMES, a split junction of the same size whose
// breakpoints both lie within this distance of the SV's is taken as this
// deletion with imprecise breakpoints, unless another deletion record of that
// size lies within this distance more of the junction (pairs of deletions too
// close to tell apart).
static const int SV_NEAR_JUNCTION = 150;

// Aligned part of an alignment: reference span [refStart, refEnd), the S/H
// clip at each end, and the query offset just past the last aligned base
// (in the orientation of the alignment).
struct SVAlignedSpan {
  int refStart = 0;
  int refEnd = 0;
  int leadClip = 0;
  int trailClip = 0;
  int queryEnd = 0;
};

static SVAlignedSpan getAlignedSpan(const bam1_t &aln) {
  const uint32_t *cigar = bam_get_cigar(&aln);
  int n = aln.core.n_cigar;
  SVAlignedSpan span;
  span.refStart = static_cast<int>(aln.core.pos);
  span.refEnd = static_cast<int>(bam_endpos(&aln));
  int first = 0;
  while (first < n && (bam_cigar_op(cigar[first]) == BAM_CSOFT_CLIP ||
                       bam_cigar_op(cigar[first]) == BAM_CHARD_CLIP))
    span.leadClip += bam_cigar_oplen(cigar[first++]);
  int last = n - 1;
  while (last >= first && (bam_cigar_op(cigar[last]) == BAM_CSOFT_CLIP ||
                           bam_cigar_op(cigar[last]) == BAM_CHARD_CLIP))
    span.trailClip += bam_cigar_oplen(cigar[last--]);
  span.queryEnd = span.leadClip;
  for (int i = first; i <= last; i++) {
    // M, I, =, X consume the query
    if (bam_cigar_type(bam_cigar_op(cigar[i])) & 1)
      span.queryEnd += bam_cigar_oplen(cigar[i]);
  }
  return span;
}

// One SA segment of the read: 0-based reference start and the query interval
// it aligns (in the orientation of the alignment).
struct SVSplitSegment {
  int refStart = 0;
  int queryStart = 0;
  int queryEnd = 0;
};

// SA segments on the same contig and strand as the alignment, with mapping
// quality >= mappingQuality, that align a later part of the read.
static std::vector<SVSplitSegment>
getDownstreamSegments(const bam1_t &aln, const SVAlignedSpan &span,
                      const std::string &chrName, int mappingQuality) {
  std::vector<SVSplitSegment> segments;
  uint8_t *saData = bam_aux_get(&aln, "SA");
  if (saData == NULL)
    return segments;
  const char *sa = bam_aux2Z(saData);
  if (sa == NULL)
    return segments;
  char strand = bam_is_rev(&aln) ? '-' : '+';
  // SA:Z:(rname,pos,strand,CIGAR,mapQ,NM;)+
  std::istringstream saStream(sa);
  std::string entry;
  while (std::getline(saStream, entry, ';')) {
    std::vector<std::string> fields;
    std::istringstream entryStream(entry);
    std::string field;
    while (std::getline(entryStream, field, ','))
      fields.push_back(field);
    if (fields.size() < 6 || fields[0] != chrName || fields[2].size() != 1 ||
        fields[2][0] != strand)
      continue;
    char *mapqEnd;
    long mapq = std::strtol(fields[4].c_str(), &mapqEnd, 10);
    if (mapqEnd == fields[4].c_str() || *mapqEnd != '\0' ||
        mapq < mappingQuality)
      continue;
    char *posEnd;
    long saPos = std::strtol(fields[1].c_str(), &posEnd, 10);
    if (posEnd == fields[1].c_str() || *posEnd != '\0' || saPos < 1 ||
        saPos > INT_MAX)
      continue;
    SVSplitSegment segment;
    segment.refStart = (int)saPos - 1;
    long long queryStart = 0;
    long long queryEnd = 0;
    bool leading = true;
    bool valid = !fields[3].empty();
    const char *op = fields[3].c_str();
    while (valid && *op) {
      char *opEnd;
      long len = std::strtol(op, &opEnd, 10);
      // a BAM CIGAR operation is at most 2^28 - 1 long
      if (opEnd == op || *opEnd == '\0' || len < 0 || len >= (1L << 28)) {
        valid = false;
        break;
      }
      char type = *opEnd;
      if (type == 'S' || type == 'H') {
        if (leading)
          queryStart += len;
      } else {
        if (leading)
          queryEnd = queryStart;
        leading = false;
        if (type == 'M' || type == 'I' || type == '=' || type == 'X')
          queryEnd += len;
      }
      op = opEnd + 1;
    }
    if (!valid || leading || queryEnd > INT_MAX)
      continue;
    segment.queryStart = static_cast<int>(queryStart);
    segment.queryEnd = static_cast<int>(queryEnd);
    if (segment.queryEnd > span.queryEnd)
      segments.push_back(segment);
  }
  return segments;
}

// Deleting ref[a, a + length) or ref[b, b + length) leaves the same sequence
// when the bases between the two starts repeat length bases later, i.e. the
// two are one deletion placed at different points of a repeat.
static bool sameDeletion(const std::string &ref, int a, int b, int length) {
  if (a > b)
    std::swap(a, b);
  if (a < 0 || length <= 0 || (size_t)b + (size_t)length > ref.size())
    return false;
  for (int i = a; i < b; i++) {
    if (std::toupper((unsigned char)ref[i]) !=
        std::toupper((unsigned char)ref[i + length]))
      return false;
  }
  return true;
}

// Distance from a record position to the left breakpoint of a split junction,
// which lies anywhere in [junctionStart, junctionEnd] when the two segments of
// the read share sequence (microhomology or a repeat).
static long long junctionDistance(long long position, long long junctionStart,
                                  long long junctionEnd) {
  if (position < junctionStart)
    return junctionStart - position;
  if (position > junctionEnd)
    return position - junctionEnd;
  return 0;
}

// Another deletion record of the SV VCF (deletionRecords: position, |SVLEN|;
// phased by phase or not), of the same size as a read's split junction, lies
// closer to the junction than distance: the read carries that deletion rather
// than a shifted copy of the SV being evaluated.
static bool otherDeletionFitsJunction(
    const std::vector<std::pair<int, int>> *deletionRecords, int svStart,
    long long junctionStart, long long junctionEnd, long long junctionDeletion,
    long long distance, double svThreshold) {
  if (deletionRecords == nullptr || distance <= 0)
    return false;
  long long low = std::min<long long>(
      std::max<long long>(junctionStart - distance + 1, INT_MIN), INT_MAX);
  std::vector<std::pair<int, int>>::const_iterator record = std::lower_bound(
      deletionRecords->begin(), deletionRecords->end(),
      std::make_pair((int)low, INT_MIN));
  for (; record != deletionRecords->end() &&
         (long long)(*record).first < junctionEnd + distance;
       record++) {
    // records at the SV's own position (the SV itself)
    if ((*record).first == svStart)
      continue;
    double size = (*record).second;
    if (std::fabs(size - (double)junctionDeletion) > svThreshold * size)
      continue;
    if (junctionDistance((*record).first, junctionStart, junctionEnd) < distance)
      return true;
  }
  return false;
}

// Allele of a deletion at VCF POS = svPos + 1 shown by a split alignment that
// ends, clipped, near the SV:
//   1  the read continues after the deletion: the next segment starts
//      |SVLEN| bases further on the reference
//   0  the read's own split deletion contains the SV position, so the read
//      shows the reference allele there (overlapping SVs)
//   2  a junction of the same size near the SV that cannot be told apart from
//      it: no evidence for either allele
//   3  a junction of the same size whose breakpoints both lie within
//      SV_NEAR_JUNCTION of the SV's, with no other deletion record that close:
//      most likely this deletion with imprecise breakpoints
//  -1  none of these
static int splitAlignmentDeletionAllele(
    const bam1_t &aln, const SVAlignedSpan &span, int svPos, int svLength,
    bool isDeletion, double svThreshold, int mappingQuality,
    const std::string &chrName, const std::string &ref_string,
    const std::vector<std::pair<int, int>> *deletionRecords) {
  // only the upstream part of a split read, clipped where it ends
  if (span.trailClip == 0 || span.refStart > svPos || svLength == INT_MIN)
    return -1;
  std::vector<SVSplitSegment> segments =
      getDownstreamSegments(aln, span, chrName, mappingQuality);
  int pos = svPos + 1;
  int delLength = std::abs(svLength);
  bool sameSizeNearby = false;
  bool nearJunction = false;
  for (const SVSplitSegment &segment : segments) {
    long long queryGap = (long long)segment.queryStart - span.queryEnd;
    long long overlap = std::max(0LL, -queryGap);
    long long insertion = std::max(0LL, queryGap);
    // a deletion between two copies of a long repeat gives segments that
    // share up to one copy of the read; beyond SV_MAX_QUERY_OVERLAP accept
    // that only for a deletion that is longer than the shared part
    if (overlap > SV_MAX_QUERY_OVERLAP && (!isDeletion || overlap >= delLength))
      continue;
    long long junctionDeletion = (long long)segment.refStart - span.refEnd + overlap;
    long long leftOffset = (long long)span.refEnd - pos;
    long long rightOffset = (long long)segment.refStart - ((long long)pos + delLength);
    bool sameSize = isDeletion && junctionDeletion > 0 &&
                    insertion <= std::max(SV_JUNCTION_TOLERANCE,
                                          (int)(svThreshold * delLength)) &&
                    (double)std::llabs(junctionDeletion - delLength) <=
                        svThreshold * delLength;
    // another deletion record of this size lies closer to the junction
    long long junctionStart = span.refEnd - overlap;
    long long ownDistance = junctionDistance(pos, junctionStart, span.refEnd);
    bool otherFitsBetter = sameSize && otherDeletionFitsJunction(
        deletionRecords, pos, junctionStart, span.refEnd, junctionDeletion,
        ownDistance, svThreshold);
    // this segment forms the junction of the deletion
    if (sameSize && !otherFitsBetter && leftOffset >= -SV_JUNCTION_TOLERANCE &&
        leftOffset <= overlap + SV_JUNCTION_TOLERANCE &&
        rightOffset >= -(overlap + SV_JUNCTION_TOLERANCE) &&
        rightOffset <= SV_JUNCTION_TOLERANCE)
      return 1;
    // the same deletion placed elsewhere in a repeat (e.g. a left-aligned
    // VCF record): its first deleted base (VCF POS, or POS - 1 for callers
    // that report the first deleted base) gives the same sequence
    if (sameSize && !otherFitsBetter && junctionDeletion == delLength &&
        std::llabs(junctionStart - pos) <=
            std::max(SV_MAX_JUNCTION_SHIFT, delLength) &&
        (sameDeletion(ref_string, (int)junctionStart, pos, delLength) ||
         sameDeletion(ref_string, (int)junctionStart, pos - 1, delLength)))
      return 1;
    // a junction of the same size that matches one breakpoint, or that is
    // shifted as a whole (breakpoints in a repeat), may still be this
    // deletion: it is not evidence for either allele
    if (sameSize && std::llabs(leftOffset) <= overlap + SV_NEAR_JUNCTION &&
        std::llabs(rightOffset) <= overlap + SV_NEAR_JUNCTION &&
        std::llabs(leftOffset) <= SV_MAX_JUNCTION_SHIFT &&
        std::llabs(rightOffset) <= SV_MAX_JUNCTION_SHIFT &&
        !otherDeletionFitsJunction(deletionRecords, pos, junctionStart,
                                   span.refEnd, junctionDeletion,
                                   ownDistance + SV_NEAR_JUNCTION,
                                   svThreshold)) {
      nearJunction = true;
      continue;
    }
    long long tolerance = overlap + SV_JUNCTION_TOLERANCE;
    if (sameSize && (std::llabs(leftOffset) <= tolerance ||
                     std::llabs(rightOffset) <= tolerance ||
                     std::llabs(leftOffset - rightOffset) <= tolerance ||
                     (std::llabs(leftOffset) <= SV_MAX_JUNCTION_SHIFT &&
                      std::llabs(rightOffset) <= SV_MAX_JUNCTION_SHIFT)) &&
        !otherFitsBetter) {
      sameSizeNearby = true;
      continue;
    }
    // the read's own deletion contains the SV position
    if (queryGap <= SV_JUNCTION_TOLERANCE &&
        (long long)segment.refStart - span.refEnd >= 2 * SV_REF_FLANK &&
        span.refEnd + SV_REF_FLANK <= pos &&
        pos <= segment.refStart - SV_REF_FLANK)
      return 0;
  }
  return nearJunction ? 3 : sameSizeNearby ? 2 : -1;
}

// Index of the CIGAR operation that consumes reference position refPos, or -1
// if the aligned part of the alignment does not reach it.
static int findReferenceCigarOp(const bam1_t &aln, int refPos) {
  const uint32_t *cigar = bam_get_cigar(&aln);
  int pos = static_cast<int>(aln.core.pos);
  for (int i = 0; i < (int)aln.core.n_cigar; i++) {
    // M, D, N, =, X consume the reference
    if (!(bam_cigar_type(bam_cigar_op(cigar[i])) & 2))
      continue;
    int len = bam_cigar_oplen(cigar[i]);
    if (refPos < pos + len)
      return refPos >= pos ? i : -1;
    pos += len;
  }
  return -1;
}

void BamParser::handleSVVariant(
    const bam1_t &aln, int cigar_i, int svWindow, double svThreshold,
    int mappingQuality, const std::string &ref_string,
    std::vector<std::pair<int, int>>::iterator &currentSVIter, int &svPos,
    ReadVariant &tmpReadResult) {

  const uint32_t *cigar = bam_get_cigar(&aln);
  int aln_core_n_cigar = aln.core.n_cigar;

  int allele = 0;
  int sv_start = (*currentSVIter).first;
  int sv_length = (*currentSVIter).second;
  int sv_end = sv_start + std::abs(sv_length);
  double sv_region = sv_end - sv_start + 1;

  // look for an I/D of the SV size around the CIGAR operation that holds the
  // SV; when the SV was reached on a clip, use the operation that covers the
  // SV instead (an insertion op is kept: it is usually the SV itself)
  int center = cigar_i;
  if (bam_cigar_op(cigar[cigar_i]) == BAM_CSOFT_CLIP ||
      bam_cigar_op(cigar[cigar_i]) == BAM_CHARD_CLIP) {
    int covering = findReferenceCigarOp(aln, svPos);
    if (covering >= 0)
      center = covering;
  }
  for (int j = std::max(center - svWindow, 0);
       j < std::min(center + svWindow, aln_core_n_cigar); j++) {
    int op = bam_cigar_op(cigar[j]);
    int oplen = bam_cigar_oplen(cigar[j]);
    if (op == BAM_CINS &&
        std::abs(sv_region - oplen) / std::abs(sv_region) < svThreshold) {
      allele = 1;
      break;
    }
    if (op == BAM_CDEL &&
        std::abs(sv_region - oplen) / std::abs(sv_region) < svThreshold) {
      allele = 1;
      break;
    }
  }

  // An alignment that does not reach well past the SV on both sides (one that
  // ends or starts at the SV, often clipped there, or that reached the SV
  // only on an S/H operation) does not show the reference allele. Decide it
  // from its split alignment or from the SV caller's RNAMES, or skip it.
  bool hasCall = true;
  if (allele == 0) {
    SVAlignedSpan span = getAlignedSpan(aln);
    int pos = svPos + 1;
    bool spans =
        span.refStart <= pos - SV_REF_FLANK && span.refEnd >= pos + SV_REF_FLANK;
    // an alignment that does reach past the SV but is clipped soon after it
    // may be the upstream part of a carrier whose alignment ran on over a
    // repeat; it keeps its REF call unless split or RNAMES evidence says ALT
    bool clippedAfter =
        span.trailClip > SV_MIN_BREAK_CLIP &&
        span.refEnd - pos <= std::max<long long>(SV_MAX_JUNCTION_SHIFT,
                                                 std::llabs((long long)sv_length));
    if (!spans || clippedAfter) {
      const SVReadSupport *support = NULL;
      if (currentSVSupport != nullptr) {
        std::map<int, SVReadSupport>::const_iterator supportIter =
            currentSVSupport->find(sv_start);
        if (supportIter != currentSVSupport->end())
          support = &(*supportIter).second;
      }
      int splitAllele = -1;
      if (support != NULL)
        splitAllele = splitAlignmentDeletionAllele(
            aln, span, svPos, sv_length, support->isDeletion, svThreshold,
            mappingQuality, chrName, ref_string, currentDeletionRecords);
      if (support != NULL && support->hasReadNames) {
        // the SV caller listed the reads that carry the SV (RNAMES): only
        // those count as ALT; a junction it did not list is no evidence
        if (span.refStart <= svPos &&
            support->readNames.count(bam_get_qname(&aln)) != 0)
          allele = 1;
        else if (splitAllele == 1 || splitAllele == 2 || splitAllele == 3 ||
                 (splitAllele == -1 && !spans))
          hasCall = false;
      } else if (splitAllele == 1 || splitAllele == 3) {
        allele = 1;
      } else if (splitAllele == 2 || (splitAllele == -1 && !spans)) {
        hasCall = false;
      }
    }
  }

  if (hasCall) {
    // use quality -1 to identify SVs
    tmpReadResult.variantVec.push_back(
        Variant(svPos, allele, -1, VariantType::SV));
  }
  currentSVIter++;
  if (currentSVIter != SV_map[chrName].end()) {
    svPos = (*currentSVIter).first - 1;
  } else {
    svPos = INT_MAX;
  }
}

void BamParser::get_snp(const bam_hdr_t &bamHdr, const bam1_t &aln,
                        std::vector<ReadVariant> &readVariantVec,
                        ClipCount &clipCount, const std::string &ref_string,
                        bool isONT, int svWindow, double svThreshold,
                        int mappingQuality) {

  ReadVariant tmpReadResult;
  tmpReadResult.read_name = bam_get_qname(&aln);
  tmpReadResult.source_id = bamHdr.target_name[aln.core.tid];
  tmpReadResult.reference_start = aln.core.pos;
  tmpReadResult.is_reverse = bam_is_rev(&aln);

  // position relative to reference
  int ref_pos = aln.core.pos;

  // position relative to read
  int query_pos = 0;

  // Skip variants that are to the left of this read
  while (firstVariantIter != currentVariants->end() &&
         (*firstVariantIter).first < ref_pos)
    firstVariantIter++;

  // Skip structure variants that are to the left of this read
  while (firstSVIter != SV_map[chrName].end() && (*firstSVIter).first < ref_pos)
    firstSVIter++;

  // Skip modify that are to the left of this read
  while (firstModIter != currentMod->end() && (*firstModIter).first < ref_pos)
    firstModIter++;

  // set variant start for current alignment
  std::map<int, SnpVariant>::iterator currentVariantIter = firstVariantIter;

  // set structure variant start for current alignment
  std::vector<std::pair<int, int>>::iterator currentSVIter = firstSVIter;

  // set modify start for current alignment
  std::map<int, std::map<std::string, ModVariant>>::iterator currentModIter =
      firstModIter;

  // set cigar pointer and number of cigar
  uint32_t *cigar = bam_get_cigar(&aln);
  int aln_core_n_cigar = aln.core.n_cigar;

  // reading cigar to detect varaint on this read
  for (int i = 0; i < aln_core_n_cigar; i++) {

    // get current cigar type and cigar length
    int cigar_op = bam_cigar_op(cigar[i]);
    int cigar_oplen = bam_cigar_oplen(cigar[i]);

    // get the starting position of each variant currently.
    // Use INT_MAX as sentinel when iterator reaches end(), to avoid
    // dereferencing end() (undefined behavior) and to ensure exhausted variant
    // types never "win" the min-position comparison in the processing loop
    // below.
    int modPos = (currentModIter != currentMod->end()) ? (*currentModIter).first
                                                       : INT_MAX;
    int svPos =
        (currentSVIter != SV_map[chrName].end()) ? (*currentSVIter).first - 1
                                                 : INT_MAX;
    int variantPos = (currentVariantIter != currentVariants->end())
                         ? (*currentVariantIter).first
                         : INT_MAX;

    // Advance currentVariantIter past any variants to the left of ref_pos.
    // NOTE (Problem⑤): This pre-advance exists only for currentVariantIter.
    // currentSVIter and currentModIter do not have a corresponding pre-advance
    // because SV and MOD variants are handled differently in the processing loop:
    // SV uses a window scan that tolerates stale positions, and MOD relies on a
    // qname lookup rather than positional ordering. The asymmetry is intentional.
    while (currentVariantIter != currentVariants->end() &&
           variantPos < ref_pos) {
      currentVariantIter++;
      if (currentVariantIter != currentVariants->end()) {
        variantPos = (*currentVariantIter).first;
      }
    }

    // Processing the region covered by the current CIGAR operator
    // Determine if any variant is included in the current CIGAR operator
    while ((currentModIter != currentMod->end() &&
            modPos < ref_pos + cigar_oplen) ||
           (currentSVIter != SV_map[chrName].end() &&
            svPos < ref_pos + cigar_oplen) ||
           (currentVariantIter != currentVariants->end() &&
            variantPos < ref_pos + cigar_oplen)) {

      // modification's position is minimal (or equal, MOD takes priority)
      if ((currentVariantIter == currentVariants->end() ||
           modPos <= variantPos) &&
          (currentSVIter == SV_map[chrName].end() || modPos <= svPos) &&
          currentModIter != currentMod->end()) {

        handleModVariant(aln, variantPos, currentModIter, modPos, tmpReadResult);
      }
      // SV's position is minimal (or equal to SNP, SV takes priority over SNP)
      else if ((currentVariantIter == currentVariants->end() ||
                svPos <= variantPos) &&
               (currentModIter == currentMod->end() || svPos < modPos) &&
               currentSVIter != SV_map[chrName].end()) {
        // If this read not contain SV, it means this read is the same as
        // reference genome. default this read the same as ref

        handleSVVariant(aln, i, svWindow, svThreshold, mappingQuality,
                        ref_string, currentSVIter, svPos, tmpReadResult);
      }

      // SNP's position is minimal
      else if ((currentSVIter == SV_map[chrName].end() || variantPos < svPos) &&
               (currentModIter == currentMod->end() || variantPos < modPos) &&
               currentVariantIter != currentVariants->end()) {

        // CIGAR operators: MIDNSHP=X correspond 012345678
        // 0: alignment match (can be a sequence match or mismatch)
        // 7: sequence match
        // 8: sequence mismatch
        if (cigar_op == 0 || cigar_op == 7 || cigar_op == 8) {
          int offset = variantPos - ref_pos;
          // The position of the variant exceeds the length of the read.
          if (query_pos + offset + 1 > int(aln.core.l_qseq)) {
            return;
          }
          handleSNPVariant(aln, i, cigar_oplen, ref_pos, query_pos,
                           currentVariantIter, variantPos, tmpReadResult);
        } else
          break;
      }
    }

    // Preparing to process the next CIGAR operator.

    // CIGAR operators: MIDNSHP=X correspond 012345678
    // 0: alignment match (can be a sequence match or mismatch)
    // 7: sequence match
    // 8: sequence mismatch
    if (cigar_op == 0 || cigar_op == 7 || cigar_op == 8) {
      query_pos += cigar_oplen;
      ref_pos += cigar_oplen;
    }
    // 1: insertion to the reference
    else if (cigar_op == 1) {
      query_pos += cigar_oplen;
    } else if (cigar_op == 2) {
      // If a reference is given, check whether any SNP falls in a homopolymer
      // region spanned by this deletion and handle it specially.
      if (ref_string != "") {
        if (handleDeletionHomopolymerSNP(aln, ref_pos, query_pos, cigar_oplen,
                                         ref_string, currentVariantIter,
                                         tmpReadResult))
          return;
      }
      ref_pos += cigar_oplen;
    }
    // 3: skipped region from the reference
    else if (cigar_op == 3) {
      ref_pos += cigar_oplen;
    }
    // 4: soft clipping (clipped sequences present in SEQ)
    else if (cigar_op == 4) {
      query_pos += cigar_oplen;
      getClip(ref_pos, i, cigar_oplen, clipCount);
    }
    // 5: hard clipping (clipped sequences NOT present in SEQ)
    else if (cigar_op == 5) {
      getClip(ref_pos, i, cigar_oplen, clipCount);
    }
    // 6: padding (silent deletion from padded reference)
    else if (cigar_op == 6) {

    } else {
      std::cerr << "alignment find unsupported CIGAR operation from read: "
                << bam_get_qname(&aln) << "\n";
      exit(1);
    }
  }

  if (tmpReadResult.variantVec.size() > 0) {
    readVariantVec.push_back(tmpReadResult);
  }
}

void BamParser::getClip(int pos, int clipFrontBack, int len,
                        ClipCount &clipCount) {
  if (len > 5) {
    if (clipFrontBack == FRONT) {
      clipCount[pos][FRONT]++;
    } else {
      clipCount[pos][BACK]++;
    }
  }
}

METHParser::METHParser(PhasingParameters &in_params, SnpParser &in_snpFile,
                       SVParser &in_svFile) {
  params = &in_params;
  snpFile = &in_snpFile;
  svFile = &in_svFile;
  representativePos = -1;
  upMethPos = -1;

  chrVariant = new std::map<std::string,
                            std::map<int, std::map<std::string, ModVariant>>>;
  representativeMap = new std::map<int, int>;

  if (params->modFile.find("gz") != std::string::npos) {
    // .vcf.gz
    compressParser(params->modFile);
  } else if (params->modFile.find("vcf") != std::string::npos) {
    // .vcf
    unCompressParser(params->modFile);
  }
}

void METHParser::writeResult(PhasingResult &phasingResult,
                             const std::string &outPrefix) {
  dispatchWriteResult(params->modFile, outPrefix + "_mod.vcf",
                      phasingResult);
}

void METHParser::writeMetaHeader(const std::string &input, bool &ps_def,
                                 std::ofstream &resultVcf) {
  if (input.substr(0, 16) == "##FORMAT=<ID=PS,") {
    ps_def = true;
  }
  // Same strategy as SVParser: piggy-back PE/H1/H2 on existing ##INFO
  // declarations. modcall's output always emits RS / MR / NR under
  // ##INFO=, so this reliably places them together.
  resultVcf << input << "\n";
  if (input.substr(0, 7) == "##INFO=") {
    writePeInfoHeaders(resultVcf);
  }
}

METHParser::~METHParser() {
  delete chrVariant;
  delete representativeMap;
}

void METHParser::parserProcess(std::string &input) {
  if (input[0] == '#')
    return;

  std::istringstream iss(input);
  std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),
                                  std::istream_iterator<std::string>());

  if (fields.size() == 0) {
    return;
  }

  // trans 1-base position to 0-base
  int pos = std::stoi(fields[1]) - 1;
  std::string chr = fields[0];

  // In a series of consecutive methylation positions,
  // the first methylation position will be used as the representative after
  // merging.
  if (upMethPos + 1 != pos) {
    representativePos = pos;
  }

  int modify_start = findFieldValueStart(fields[8], fields[9], "GT");

  // homo GT
  if (fields[9][modify_start] == fields[9][modify_start + 2]) {
    return;
  }

  // conflict pos with SNP and SV
  if ((*snpFile).findSNP(chr, pos) || (*svFile).findSV(chr, pos)) {
    return;
  }

  bool is_reverse;
  // get strand
  if (fields[7].find("RS=P") != std::string::npos) {
    is_reverse = false;
  } else if (fields[7].find("RS=N") != std::string::npos) {
    is_reverse = true;
  } else {
    return;
  }

  // parse MR and NR reads
  parseReadList(fields[7], "MR=", chr, representativePos, is_reverse, true,  *chrVariant);
  parseReadList(fields[7], "NR=", chr, representativePos, is_reverse, false, *chrVariant);

  // Record the positions corresponding to the current representative position
  (*representativeMap)[pos] = representativePos;
  upMethPos = pos;
}

void METHParser::writeDataLine(const std::string &input,
                                   std::ofstream &resultVcf,
                                   PhasingResult &phasingResult) {
  std::istringstream iss(input);
  std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),
                                  std::istream_iterator<std::string>());

  if (fields.size() == 0)
    return;

  int pos = std::stoi(fields[1]);
  // use current position to find representative position
  int posIdx = (*representativeMap)[pos - 1];

  std::string key = fields[0] + "_" + std::to_string(posIdx);

  PhasingResult::iterator psElementIter = phasingResult.find(key);

  // PS flag already exist, erase PS info
  if (fields[8].find("PS") != std::string::npos) {

    std::string format_before_ps_erase = fields[8];

    // erase PS flag
    int ps_pos = fields[8].find("PS");
    if (fields[8].find(":", ps_pos + 1) != std::string::npos) {
      fields[8].erase(ps_pos, 3);
    } else {
      fields[8].erase(ps_pos - 1, 3);
    }

    // find PS value start
    int ps_start = findFieldValueStart(format_before_ps_erase, fields[9], "PS");

    // erase PS value
    if (fields[9].find(":", ps_start + 1) != std::string::npos) {
      int ps_end_pos = fields[9].find(":", ps_start + 1);
      fields[9].erase(ps_start, ps_end_pos - ps_start + 1);
    } else {
      fields[9].erase(ps_start - 1, fields[9].length() - ps_start + 1);
    }
  }
  // reset GT flag
  if (fields[8].find("GT") != std::string::npos) {
    int modify_start = findFieldValueStart(fields[8], fields[9], "GT");
    if (fields[9][modify_start + 1] == '|') {
      // direct modify GT value
      if (fields[9][modify_start] > fields[9][modify_start + 2]) {
        fields[9][modify_start + 1] = fields[9][modify_start];
        fields[9][modify_start] = fields[9][modify_start + 2];
        fields[9][modify_start + 2] = fields[9][modify_start + 1];
      }
      fields[9][modify_start + 1] = '/';
    }
  }

  // Check if the variant is extracted from this VCF
  auto posIter = (*chrVariant)[fields[0]].find(posIdx);

  // this pos is phase and exist in map
  if (psElementIter != phasingResult.end() &&
      posIter != (*chrVariant)[fields[0]].end()) {

    // add PS flag and value
    fields[8] = fields[8] + ":PS";
    fields[9] =
        fields[9] + ":" + std::to_string((*psElementIter).second.block);

    int modify_start = findFieldValueStart(fields[8], fields[9], "GT");
    // direct modify GT value
    fields[9][modify_start] = (*psElementIter).second.RAstatus[0];
    fields[9][modify_start + 1] = '|';
    fields[9][modify_start + 2] = (*psElementIter).second.RAstatus[2];

    // Attach PE / H1 / H2 to INFO for phased positions only.
    // Rule: if INFO was ".", replace it; otherwise append with ';'.
    appendPeHapInfoToInfoField(fields[7], (*psElementIter).second);
  }
  // this pos has not been phased
  else {
    // drop stale PE/H1/H2 left over from a previously phased input
    removePeHapFromInfoField(fields[7]);
    // add PS flag and value
    fields[8] = fields[8] + ":PS";
    fields[9] = fields[9] + ":.";
  }
  for (std::vector<std::string>::iterator fieldIter = fields.begin();
       fieldIter != fields.end(); ++fieldIter) {
    if (fieldIter != fields.begin())
      resultVcf << "\t";
    resultVcf << (*fieldIter);
  }
  resultVcf << "\n";
}

std::map<int, std::map<std::string, ModVariant>>
METHParser::getVariants(std::string chrName) {
  std::map<int, std::map<std::string, ModVariant>> targetVariants;
  std::map<std::string,
           std::map<int, std::map<std::string, ModVariant>>>::iterator chrIter =
      chrVariant->find(chrName);

  if (chrIter != chrVariant->end())
    targetVariants = (*chrIter).second;

  return targetVariants;
}
