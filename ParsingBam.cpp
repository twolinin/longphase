#include "ParsingBam.h"

#include <cstdio>
#include <cmath>
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

    // read file
    std::string chr_info(faidx_fetch_seq(fai, (*iter).c_str(), 0,
                                         last_pos.at(index) + 5, &ref_len));
    if (ref_len == 0) {
      std::cout << "nothing in reference file \n";
    }

    // update map
    chrString[(*iter)] = chr_info;
  }
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
                                      PhasingResult phasingResult) {
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
                                        PhasingResult phasingResult) {
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
                                            PhasingResult phasingResult) {
  if (inputFile.find("gz") != std::string::npos)
    compressInput(inputFile, outputFile, phasingResult);
  else if (inputFile.find("vcf") != std::string::npos)
    unCompressInput(inputFile, outputFile, phasingResult);
}

void BaseVairantParser::writeLine(std::string &input, bool &ps_def,
                                  std::ofstream &resultVcf,
                                  PhasingResult &phasingResult) {
  if (input.substr(0, 2) == "##")
    writeMetaHeader(input, ps_def, resultVcf);
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
    resultVcf << "##longphaseVersion=" << params->version << "\n";
    resultVcf << "##commandline=\"" << params->command << "\"\n";
    commandLine = true;
  }
  resultVcf << input << "\n";
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
  // set all sample string
  std::string allSmples = "-";
  // limit the VCF data to the sample name passed in
  int is_file = bcf_hdr_set_samples(hdr, allSmples.c_str(), 0);
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

SnpParser::~SnpParser() { delete chrVariant; }

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

void SnpParser::writeResult(PhasingResult phasingResult) {
  dispatchWriteResult(params->snpFile, params->resultPrefix + ".vcf",
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

  // this pos is phase
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
  }
  // this pos has not been phased
  else {
    // add PS flag and value
    fields[8] = fields[8] + ":PS";
    fields[9] = fields[9] + ":.";
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

void SVParser::writeResult(PhasingResult phasingResult) {
  dispatchWriteResult(params->svFile, params->resultPrefix + "_SV.vcf",
                      phasingResult);
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
  }
  // this pos has not been phased
  else {
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
            params.isONT, params.svWindow, params.svThreshold);
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

void BamParser::handleSVVariant(
    const bam1_t &aln, int cigar_i, int svWindow, double svThreshold,
    std::vector<std::pair<int, int>>::iterator &currentSVIter, int &svPos,
    ReadVariant &tmpReadResult) {

  const uint32_t *cigar = bam_get_cigar(&aln);
  int aln_core_n_cigar = aln.core.n_cigar;

  int allele = 0;
  int sv_start = (*currentSVIter).first;
  int sv_length = (*currentSVIter).second;
  int sv_end = sv_start + std::abs(sv_length);
  double sv_region = sv_end - sv_start + 1;

  for (int j = std::max(cigar_i - svWindow, 0);
       j < std::min(cigar_i + svWindow, aln_core_n_cigar); j++) {
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

  // use quality -1 to identify SVs
  tmpReadResult.variantVec.push_back(
      Variant(svPos, allele, -1, VariantType::SV));
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
                        bool isONT, int svWindow, double svThreshold) {

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

        handleSVVariant(aln, i, svWindow, svThreshold, currentSVIter, svPos, tmpReadResult);
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

void METHParser::writeResult(PhasingResult phasingResult) {
  dispatchWriteResult(params->modFile, params->resultPrefix + "_mod.vcf",
                      phasingResult);
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
  }
  // this pos has not been phased
  else {
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
