#include "ModCallParsingBam.h"
#include "PhasingGraph.h"
#include "Util.h"

#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <typeinfo>

MethFastaParser::MethFastaParser(std::string fastaFile,
                                 std::vector<ReferenceChromosome> &chrInfo) {
  if (fastaFile == "") {
    return;
  }

  faidx_t *fai = fai_load(fastaFile.c_str());
  int fai_nseq = faidx_nseq(fai);
  const char *seqname;
  int seqlen = 0;
  for (int i = 0; i < fai_nseq; i++) {
    int ref_len = 0;
    seqname = faidx_iseq(fai, i);
    seqlen = faidx_seq_len(fai, seqname);
    char *seq = faidx_fetch_seq(fai, seqname, 0, seqlen + 1, &ref_len);
    chrInfo.emplace_back(seqname, seq, seqlen);
    free(seq);
  }
  fai_destroy(fai);
}

MethFastaParser::~MethFastaParser() {}

MethBamParser::MethBamParser(std::string inputChrName,
                             ModCallParameters &in_params,
                             MethSnpParser &snpMap, std::string &ref_string)
    : chrName(inputChrName) {
  params = &in_params;
  refString = &ref_string;
  refstartpos = 0;
  chrMethMap = new std::map<int, MethPosInfo>;
  readStartEndMap = new std::map<int, std::pair<int, int>>;

  currentVariants = new std::map<int, SnpVariant>;

  if (snpMap.hasValidSnpData()) {
    (*currentVariants) = snpMap.getVariants_markindel(chrName, ref_string);
  }

  firstVariantIter = currentVariants->begin();
}

MethBamParser::~MethBamParser() {
  delete chrMethMap;
  delete readStartEndMap;
  delete currentVariants;
}

bool MethBamParser::isValidAlignment(const bam1_t *aln) const {
  // Intentionally different from phase/haplotag filters: fixed qual>=1;
  // supplementary (0x800) alignments are excluded.
  int flag = aln->core.flag;
  return aln->core.qual >= 1
         && (flag & 0x4) == 0   // not unmapped
         && (flag & 0x100) == 0 // not secondary
         && (flag & 0x400) == 0 // not duplicate
         && (flag & 0x800) == 0;
}

bool MethBamParser::extractSnpAllele(
    std::map<int, SnpVariant>::iterator &currentVariantIter, int ref_pos,
    int length, int querypos, int cigaridx, int aln_core_n_cigar,
    const uint32_t *cigar, const bam1_t &aln, ReadVariant *tmpReadResult) {
  while (currentVariantIter != currentVariants->end() &&
         currentVariantIter->first < ref_pos + length) {
    int variantPos = currentVariantIter->first;
    if (variantPos >= ref_pos) {
      int refAlleleLen = currentVariantIter->second.Ref.length();
      int altAlleleLen = currentVariantIter->second.Alt.length();
      int offset = variantPos - ref_pos;
      int base_q = 0;
      int allele = -1;

      if (querypos + offset + 1 > int(aln.core.l_qseq)) {
        return false;
      }

      if (refAlleleLen == 1 && altAlleleLen == 1) {
        char base =
            seq_nt16_str[bam_seqi(bam_get_seq(&aln), querypos + offset)];
        if (base == currentVariantIter->second.Ref[0])
          allele = 0;
        else if (base == currentVariantIter->second.Alt[0])
          allele = 1;

        base_q = bam_get_qual(&aln)[querypos + offset];
      }

      if (refAlleleLen == 1 && altAlleleLen != 1 &&
          cigaridx + 1 < aln_core_n_cigar) {
        if (ref_pos + length - 1 == variantPos &&
            bam_cigar_op(cigar[cigaridx + 1]) == 1) {
          allele = 1;
        } else {
          allele = 0;
        }
        base_q = -4;

        if (currentVariantIter->second.is_danger) {
          base_q = -5;
        }
      }

      if (refAlleleLen != 1 && altAlleleLen == 1 &&
          cigaridx + 1 < aln_core_n_cigar) {
        if (ref_pos + length - 1 == variantPos &&
            bam_cigar_op(cigar[cigaridx + 1]) == 2) {
          allele = 1;
        } else {
          allele = 0;
        }
        base_q = -4;

        if (currentVariantIter->second.is_danger) {
          base_q = -5;
        }
      }

      if (allele != -1) {
        tmpReadResult->variantVec.emplace_back(variantPos, allele, base_q,
                                               VariantType::SNP);
        (*chrMethMap)[variantPos].variantType = VariantType::SNP;
      }
    }
    currentVariantIter++;
  }
  return true;
}

void MethBamParser::classifyMethylationBase(
    int &n, int &pos, hts_base_mod *mods, int querypos, int length,
    int refpos, const bam1_t &aln, ReadVariant *tmpReadResult,
    hts_base_mod_state *mod_state) {
  while (true) {
    if (pos > (querypos + length)) {
      break;
    }
    if (n <= 0) {
      break;
    }

    int methrpos;
    if (bam_is_rev(&aln)) {
      methrpos = pos - querypos + refpos - 1;
    } else {
      methrpos = pos - querypos + refpos;
    }

    if (int((*refString).length()) < methrpos) {
      break;
    }

    for (int j = 0; j < n && j < 5; j++) {
      auto it = chrMethMap->find(methrpos);
      if (mods[j].modified_base == 109 && pos <= (querypos + length) &&
          (it == chrMethMap->end() ||
           it->second.variantType == VariantType::MOD)) {
        if (mods[j].qual >= params->modThreshold * 255) {
          (*chrMethMap)[methrpos].methreadcnt++;
          (*chrMethMap)[methrpos].variantType = VariantType::MOD;
          (*chrMethMap)[methrpos].strand = (bam_is_rev(&aln) ? 1 : 0);
          (*chrMethMap)[methrpos].modReadVec.emplace_back(bam_get_qname(&aln));
          tmpReadResult->variantVec.emplace_back(methrpos, 0, 60,
                                                 VariantType::MOD);
        } else if (mods[j].qual <= params->unModThreshold * 255) {
          (*chrMethMap)[methrpos].canonreadcnt++;
          (*chrMethMap)[methrpos].variantType = VariantType::MOD;
          (*chrMethMap)[methrpos].nonModReadVec.emplace_back(
              bam_get_qname(&aln));
          tmpReadResult->variantVec.emplace_back(methrpos, 1, 60,
                                                 VariantType::MOD);
        } else {
          (*chrMethMap)[methrpos].noisereadcnt++;
          (*chrMethMap)[methrpos].variantType = VariantType::MOD;
        }
      }
    }
    n = bam_next_basemod(&aln, mod_state, mods, 5, &pos);
  }
}

bool MethBamParser::handleDeletionInHomopolymer(
    std::map<int, SnpVariant>::iterator &currentVariantIter, int ref_pos,
    int length, int querypos, const bam1_t &aln,
    hts_base_mod_state *mod_state, ReadVariant *tmpReadResult) {
  int del_len = length;
  if (ref_pos + del_len + 1 == currentVariantIter->first) {
    return true;
  }

  if (currentVariantIter->first < ref_pos ||
      currentVariantIter->first >= ref_pos + del_len) {
    return true;
  }

  if (homopolymerLength(currentVariantIter->first, *refString) < 3) {
    return true;
  }

  int refAlleleLen = currentVariantIter->second.Ref.length();
  int altAlleleLen = currentVariantIter->second.Alt.length();
  int base_q = 0;

  if (querypos + 1 > aln.core.l_qseq) {
    hts_base_mod_state_free(mod_state);
    return false;
  }

  int allele = -1;
  if (refAlleleLen == 1 && altAlleleLen == 1) {
    char base = seq_nt16_str[bam_seqi(bam_get_seq(&aln), querypos)];
    if (base == currentVariantIter->second.Ref[0]) {
      allele = 0;
    } else if (base == currentVariantIter->second.Alt[0]) {
      allele = 1;
    }
    base_q = bam_get_qual(&aln)[querypos];
  } else if (refAlleleLen != 1 && altAlleleLen == 1) {
    allele = 1;
    base_q = -4;
  }

  if (allele != -1) {
    tmpReadResult->variantVec.emplace_back(currentVariantIter->first, allele,
                                           base_q, VariantType::SNP);
    (*chrMethMap)[currentVariantIter->first].variantType = VariantType::SNP;
    currentVariantIter++;
  }

  return true;
}

void MethBamParser::updateDepthBoundary(int refstart, int refpos,
                                        bool is_reverse) {
  int refend = (is_reverse ? refpos : refpos + 1);

  if (is_reverse) {
    (*readStartEndMap)[refstart + 1].second += 1;
    (*readStartEndMap)[refend].second -= 1;
  } else {
    (*readStartEndMap)[refstart + 1].first += 1;
    (*readStartEndMap)[refend].first -= 1;
  }
}

void MethBamParser::detectMeth(std::string chrName, int lastSNPPos,
                               htsThreadPool &threadPool,
                               std::vector<ReadVariant> &readVariantVec) {
  // record SNP start iter
  std::map<int, SnpVariant>::iterator tmpFirstVariantIter = firstVariantIter;
  for (auto bamFile : params->bamFileVec) {
    firstVariantIter = tmpFirstVariantIter;
    samFile *fp_in = hts_open(bamFile.c_str(), "r");
    hts_set_fai_filename(fp_in, params->fastaFile.c_str());
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    bam1_t *aln = bam_init1();
    hts_idx_t *idx = sam_index_load(fp_in, bamFile.c_str());

    if (idx == nullptr) {
      std::cout << "ERROR: Cannot open index for bam file\n";
      exit(1);
    }

    std::string range = chrName + ":1-" + std::to_string(lastSNPPos);
    hts_itr_t *iter = sam_itr_querys(idx, bamHdr, range.c_str());

    hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &threadPool);

    while (sam_itr_next(fp_in, iter, aln) >= 0) {
      if (!isValidAlignment(aln)) {
        continue;
      }
      parse_CIGAR(*bamHdr, *aln, readVariantVec);
    }

    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
    sam_close(fp_in);
  }
}

void MethBamParser::parse_CIGAR(const bam_hdr_t &bamHdr, const bam1_t &aln,
                                std::vector<ReadVariant> &readVariantVec) {
  // see detail https://github.com/samtools/htslib/blob/develop/sam_mods.c

  // hts_base_mod_state can parse MM, ML and MN tags
  hts_base_mod_state *mod_state = hts_base_mod_state_alloc();
  if (bam_parse_basemod(&aln, mod_state) < 0) {
    hts_base_mod_state_free(mod_state);
    return;
  }
  ReadVariant *tmpReadResult = new ReadVariant();
  (*tmpReadResult).read_name = bam_get_qname(&aln);
  (*tmpReadResult).source_id = bamHdr.target_name[aln.core.tid];
  (*tmpReadResult).reference_start = aln.core.pos;
  (*tmpReadResult).is_reverse = bam_is_rev(&aln);

  int refstart = aln.core.pos;
  // forward strand is checked from head to tail
  // reverse strand is checked from tail to head
  int refpos = (bam_is_rev(&aln) ? refstart + 1 : refstart);
  int ref_pos = aln.core.pos;
  // auto qmethiter = (align.is_reverse ? align.queryMethVec.end()-1 :
  // align.queryMethVec.begin() );
  int querypos = 0;
  int aln_core_n_cigar = aln.core.n_cigar;
  uint32_t *cigar = bam_get_cigar(&aln);

  hts_base_mod mods[5];
  int pos;
  // bam_next_basemod can iterate over this cached state.
  int n = bam_next_basemod(&aln, mod_state, mods, 5, &pos);
  if (n <= 0) {
    hts_base_mod_state_free(mod_state);
    delete tmpReadResult;
    return;
  }

  while (firstVariantIter != currentVariants->end() &&
         (*firstVariantIter).first < ref_pos) {
    firstVariantIter++;
  }
  // set variant start for current alignment
  std::map<int, SnpVariant>::iterator currentVariantIter = firstVariantIter;

  // Parse CIGAR
  for (int cigaridx = 0; cigaridx < int(aln.core.n_cigar); cigaridx++) {

    int cigar_op = bam_cigar_op(cigar[cigaridx]);
    int length = bam_cigar_oplen(cigar[cigaridx]);
    if (cigar_op == 0 || cigar_op == 7 || cigar_op == 8) {
      if (!extractSnpAllele(currentVariantIter, ref_pos, length, querypos,
                            cigaridx, aln_core_n_cigar, cigar, aln,
                            tmpReadResult)) {
        return;
      }
    }
    // else break;
    //}

    // CIGAR operators: MIDNSHP=X correspond 012345678
    // 0: alignment match (can be a sequence match or mismatch)
    // 7: sequence match
    // 8: sequence mismatch

    if (cigar_op == 0 || cigar_op == 7 || cigar_op == 8) {
      classifyMethylationBase(n, pos, mods, querypos, length, refpos, aln,
                              tmpReadResult, mod_state);
      querypos += length;
      refpos += length;
      ref_pos += length;
    }
    // 1: insertion to the reference
    else if (cigar_op == 1) {
      while (n > 0 && pos <= (querypos + length)) {
        n = bam_next_basemod(&aln, mod_state, mods, 5, &pos);
      }
      querypos += length;
    }
    // 2: deletion from the reference
    else if (cigar_op == 2) {
      if (*refString != "") {
        if (currentVariantIter == currentVariants->end()) {
          refpos += length;
          ref_pos += length;
          continue;
        }
        if (!handleDeletionInHomopolymer(currentVariantIter, ref_pos, length,
                                         querypos, aln, mod_state,
                                         tmpReadResult)) {
          delete tmpReadResult;
          return;
        }
      }
      refpos += length;
      ref_pos += length;
    }
    // 3: skipped region from the reference
    else if (cigar_op == 3) {
      refpos += length;
      ref_pos += length;
    }
    // 4: soft clipping (clipped sequences present in SEQ)
    else if (cigar_op == 4) {
      while (n > 0 && pos <= (querypos + length)) {
        n = bam_next_basemod(&aln, mod_state, mods, 5, &pos);
      }
      querypos += length;
    }
    // 5: hard clipping (clipped sequences NOT present in SEQ)
    // 6: padding (silent deletion from padded reference)
    else if (cigar_op == 5 || cigar_op == 6) {
      // do nothing
      refpos += length;
    }
  }

  hts_base_mod_state_free(mod_state);

  updateDepthBoundary(refstart, refpos, bam_is_rev(&aln));

  if ((*tmpReadResult).variantVec.size() > 0) {
    (*tmpReadResult).sort();
    readVariantVec.emplace_back((*tmpReadResult));
  }
  delete tmpReadResult;
}

std::string MethBamParser::buildInfoStr(const MethPosInfo &info) const {
  std::string infostr;
  if (info.modReadVec.size() > 0) {
    infostr += "MR=";
    for (const auto &readName : info.modReadVec) {
      infostr += readName + ",";
    }
    infostr.back() = ';';
  }

  if (info.nonModReadVec.size() > 0) {
    infostr += "NR=";
    for (const auto &readName : info.nonModReadVec) {
      infostr += readName + ",";
    }
    infostr.back() = ';';
  }

  return infostr;
}

std::string MethBamParser::formatMethVcfLine(
    int pos, const MethPosInfo &info, const std::string &chrSequence,
    int chrLen) const {
  if (chrLen < pos)
    return "";

  std::string ref = chrSequence.substr(pos, 1);
  if (ref != "A" && ref != "T" && ref != "C" && ref != "G" &&
      ref != "a" && ref != "t" && ref != "c" && ref != "g")
    return "";

  std::string strandinfo;
  if (info.strand == 1)
    strandinfo = "RS=N;";
  else if (info.strand == 0)
    strandinfo = "RS=P;";
  else
    return "";

  std::string infostr = buildInfoStr(info);
  std::string samplestr = info.heterstatus + ":" +
                          std::to_string(info.methreadcnt) + ":" +
                          std::to_string(info.canonreadcnt) + ":" +
                          std::to_string(info.depth);

  return chrName + "\t" + std::to_string(pos + 1) + "\t" + "." + "\t" + ref +
         "\t" + "N" + "\t" + "." + "\t" + "PASS" + "\t" + strandinfo +
         infostr + "\t" + "GT:MD:UD:DP" + "\t" + samplestr + "\n";
}

void MethBamParser::exportResult(std::string chrName, std::string chrSquence,
                                 int chrLen, std::vector<int> &passPosition,
                                 std::ostringstream &modCallResult) {

  if (params->outputAllMod) {
    for (const auto &pair : *chrMethMap) {
      std::string line =
          formatMethVcfLine(pair.first, pair.second, chrSquence, chrLen);
      if (line.empty())
        return;
      modCallResult << line;
    }
  } else {
    // used to track processed positions, avoid duplicate output
    std::set<int> processedPositions;
    // passPosition is already sorted, no need to sort again
    for (const auto &pos : passPosition) {
      // if the position is already processed, skip
      if (processedPositions.count(pos) > 0) {
        continue;
      }

      // process position pos
      auto posinfoIter = chrMethMap->find(pos);
      if (posinfoIter != chrMethMap->end()) {
        // output all modified positions or only heterozygous positions
        if (params->outputAllMod || posinfoIter->second.heterstatus == "0/1") {
          std::string line =
              formatMethVcfLine(pos, posinfoIter->second, chrSquence, chrLen);
          if (line.empty())
            continue;
          modCallResult << line;
        }
      }
      processedPositions.insert(pos);

      // process position pos+1
      int nextPos = pos + 1;
      auto nextPosinfoIter = chrMethMap->find(nextPos);
      if (nextPosinfoIter != chrMethMap->end() &&
          processedPositions.count(nextPos) == 0) {
        // output all modified positions or only heterozygous positions
        if (params->outputAllMod ||
            nextPosinfoIter->second.heterstatus == "0/1") {
          std::string line = formatMethVcfLine(
              nextPos, nextPosinfoIter->second, chrSquence, chrLen);
          if (line.empty())
            continue;
          modCallResult << line;
        }
        processedPositions.insert(nextPos);
      }
    }
  }
}

void writeResultVCF(
    ModCallParameters &params, std::vector<ReferenceChromosome> &chrInfo,
    std::map<std::string, std::ostringstream> &chrModCallResult) {

  std::ofstream modCallResultVcf(params.resultPrefix + ".vcf",
                                 std::ios_base::out | std::ios_base::trunc);
  if (!modCallResultVcf.is_open()) {
    std::cerr << "Fail to open output file :\n";
  } else {
    // set vcf header
    modCallResultVcf << "##fileformat=VCFv4.2\n";
    modCallResultVcf
        << "##INFO=<ID=RS,Number=.,Type=String,Description=\"Read Strand\">\n";
    modCallResultVcf << "##INFO=<ID=MR,Number=.,Type=String,Description=\"Read "
                        "Name of Modified position\">\n";
    modCallResultVcf << "##INFO=<ID=NR,Number=.,Type=String,Description=\"Read "
                        "Name of nonModified position\">\n";
    modCallResultVcf
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    modCallResultVcf << "##FORMAT=<ID=MD,Number=1,Type=Integer,Description="
                        "\"Modified Depth\">\n";
    modCallResultVcf << "##FORMAT=<ID=UD,Number=1,Type=Integer,Description="
                        "\"Unmodified Depth\">\n";
    modCallResultVcf << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description="
                        "\"Read Depth\">\n";
    for (const auto &chrIter : chrInfo) {
      modCallResultVcf << "##contig=<ID=" << chrIter.name
                       << ",length=" << chrIter.length << ">\n";
    }
    modCallResultVcf << "##longphaseVersion=" << params.version << "\n";
    modCallResultVcf << "##commandline=\"" << params.command << "\"\n";
    modCallResultVcf
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

    // set vcf body
    for (const auto &chrIter : chrInfo) {
      modCallResultVcf << chrModCallResult[chrIter.name].str();
    }
  }
}

void MethBamParser::computeSinglePositionGenotype() {
  for (std::map<int, MethPosInfo>::iterator chrmethmapIter =
           chrMethMap->begin();
       chrmethmapIter != chrMethMap->end(); chrmethmapIter++) {

    float methcnt = (*chrmethmapIter).second.methreadcnt;
    float nonmethcnt = (*chrmethmapIter).second.canonreadcnt;
    float depth = (*chrmethmapIter).second.depth;
    float noisereadcnt = depth - methcnt - nonmethcnt;

    if (methcnt < 0 || nonmethcnt < 0) {
      continue;
    }
    if (std::max(methcnt, nonmethcnt) == 0) {
      continue;
    }

    float heterRatio =
        std::min(methcnt, nonmethcnt) / std::max(methcnt, nonmethcnt);
    float noiseRatio = noisereadcnt / depth;

    if (heterRatio >= params->heterRatio && noiseRatio <= params->noiseRatio) {
      (*chrmethmapIter).second.heterstatus = "0/1";
    } else if (methcnt >= nonmethcnt) {
      (*chrmethmapIter).second.heterstatus = "1/1";
    } else {
      (*chrmethmapIter).second.heterstatus = "0/0";
    }
  }
}

std::set<int> MethBamParser::applyCpGStrandPairGenotype() {
  std::set<int> positionPairs;

  for (auto iter = chrMethMap->begin(); iter != chrMethMap->end(); iter++) {
    if (iter->second.strand == 0 &&
        iter->second.variantType == VariantType::MOD) {
      int currPos = iter->first;
      int nextPos = currPos + 1;

      auto nextIter = chrMethMap->find(nextPos);
      if (nextIter != chrMethMap->end() && nextIter->second.strand == 1 &&
          nextIter->second.variantType == VariantType::MOD) {
        float totalMethCnt =
            iter->second.methreadcnt + nextIter->second.methreadcnt;
        float totalNonMethCnt =
            iter->second.canonreadcnt + nextIter->second.canonreadcnt;
        float totalDepth = iter->second.depth + nextIter->second.depth;
        float totalNoiseCnt = totalDepth - totalMethCnt - totalNonMethCnt;

        if (std::max(totalMethCnt, totalNonMethCnt) == 0) {
          continue;
        }

        float combinedHeterRatio = std::min(totalMethCnt, totalNonMethCnt) /
                                   std::max(totalMethCnt, totalNonMethCnt);
        float combinedNoiseRatio = totalNoiseCnt / totalDepth;

        std::string combinedStatus;
        if (combinedHeterRatio >= params->heterRatio &&
            combinedNoiseRatio <= params->noiseRatio) {
          combinedStatus = "0/1";
          positionPairs.insert(currPos);
        } else if (totalMethCnt >= totalNonMethCnt) {
          combinedStatus = "1/1";
        } else {
          combinedStatus = "0/0";
        }
        iter->second.heterstatus = combinedStatus;
        nextIter->second.heterstatus = combinedStatus;
      }
    }
  }

  return positionPairs;
}

void MethBamParser::filterModReadVariants(
    const std::set<int> &positionPairs,
    std::vector<ReadVariant> &readVariantVec,
    std::vector<ReadVariant> &modReadVariantVec) {
  for (auto &read : readVariantVec) {
    ReadVariant newRead;
    newRead.read_name = read.read_name;
    newRead.is_reverse = read.is_reverse;

    for (auto &variant : read.variantVec) {
      int pos = variant.position;
      if (variant.type == VariantType::MOD) {
        auto methPosIter = chrMethMap->find(pos);
        if (methPosIter == chrMethMap->end()) {
          continue;
        }
        if (methPosIter->second.strand == 0) {
          if (positionPairs.count(pos)) {
            newRead.variantVec.emplace_back(pos, variant.allele,
                                            variant.quality, VariantType::MOD);
          }
        } else if (methPosIter->second.strand == 1) {
          // normalize reverse-strand position to forward-strand CpG coordinate
          if (positionPairs.count(pos - 1)) {
            newRead.variantVec.emplace_back(pos - 1, variant.allele,
                                            variant.quality, VariantType::MOD);
          }
        }
      } else if (variant.type == VariantType::SNP) {
        newRead.variantVec.emplace_back(variant);
      }
    }
    if (!newRead.variantVec.empty()) {
      modReadVariantVec.push_back(newRead);
    }
    newRead.variantVec.clear();
    newRead.variantVec.shrink_to_fit();
  }
}

void MethBamParser::judgeMethGenotype(
    std::string chrName, std::vector<ReadVariant> &readVariantVec,
    std::vector<ReadVariant> &modReadVariantVec) {
  computeSinglePositionGenotype();
  std::set<int> positionPairs = applyCpGStrandPairGenotype();
  filterModReadVariants(positionPairs, readVariantVec, modReadVariantVec);
}

void MethBamParser::dumpMethMap(const std::string &chrName,
                                std::ofstream &out) const {
  for (const auto &entry : *chrMethMap) {
    out << chrName << "\t" << entry.first << "\t"
        << entry.second.methreadcnt << "\t" << entry.second.canonreadcnt
        << "\t" << entry.second.noisereadcnt << "\t"
        << static_cast<int>(entry.second.variantType) << "\t"
        << entry.second.strand << "\n";
  }
}

void MethBamParser::dumpReadStartEnd(const std::string &chrName,
                                     std::ofstream &out) const {
  for (const auto &entry : *readStartEndMap) {
    out << chrName << "\t" << entry.first << "\t" << entry.second.first
        << "\t" << entry.second.second << "\n";
  }
}

void MethBamParser::dumpMethGenoMap(const std::string &chrName,
                                    std::ofstream &out) const {
  for (const auto &entry : *chrMethMap) {
    out << chrName
        << "\t" << entry.first
        << "\t" << entry.second.heterstatus
        << "\n";
  }
}

void MethBamParser::calculateDepth() {
  std::map<int, MethPosInfo>::iterator methIter = chrMethMap->begin();
  std::pair<int, int> currdepth = std::make_pair(0, 0);

  for (auto ReadIter = readStartEndMap->begin();
       ReadIter != readStartEndMap->end(); ReadIter++) {

    std::map<int, std::pair<int, int>>::iterator nextReadIter =
        std::next(ReadIter, 1);
    if (methIter == chrMethMap->end()) {
      break;
    }
    if (nextReadIter == readStartEndMap->end()) {
      break;
    }
    currdepth.first += (*ReadIter).second.first;
    currdepth.second += (*ReadIter).second.second;

    while (methIter != chrMethMap->end() &&
           (*methIter).first >= (*ReadIter).first &&
           (*methIter).first < (*nextReadIter).first) {

      // strand +
      if ((*methIter).second.strand == 0) {
        (*methIter).second.depth = currdepth.first;
      }
      // strand -
      else if ((*methIter).second.strand == 1) {
        (*methIter).second.depth = currdepth.second;
      }

      methIter++;
    }
  }

  readStartEndMap->clear();
}

MethylationGraph::MethylationGraph(ModCallParameters &in_params) {
  params = &in_params;
  nodeInfo = new std::map<int, std::map<std::string, VariantType>>;
  edgeList = new std::map<int, VariantEdge *>;
  forwardModNode = new std::map<int, ReadBaseMap *>;
  reverseModNode = new std::map<int, ReadBaseMap *>;
}

MethylationGraph::~MethylationGraph() {}

void MethylationGraph::addEdge(std::vector<ReadVariant> &in_readVariant,
                               std::string chrName) {
  readVariant = &in_readVariant;
  int readCount = 0;
  int vectorSize = 0;
  // iter all read
  for (std::vector<ReadVariant>::iterator readIter = in_readVariant.begin();
       readIter != in_readVariant.end(); readIter++) {
    ReadVariant tmpRead;
    vectorSize += (*readIter).variantVec.size();
    readCount++;

    for (auto variant : (*readIter).variantVec) {
      auto nodeIter = nodeInfo->find(variant.position);

      if (nodeIter == nodeInfo->end()) {
        (*nodeInfo)[variant.position] = std::map<std::string, VariantType>();
      }
      (*nodeInfo)[variant.position][(*readIter).read_name] = variant.type;
    }

    for (auto variant1Iter = (*readIter).variantVec.begin();
         variant1Iter != (*readIter).variantVec.end(); ++variant1Iter) {
      auto variant2Iter = std::next(variant1Iter);
      int searchCount = 0;
      while (variant2Iter != (*readIter).variantVec.end() && searchCount < 50) {
        if (!(variant1Iter->type == VariantType::SNP &&
              variant2Iter->type == VariantType::SNP)) {

          int pos = variant1Iter->position;
          auto edgeIter = edgeList->find(pos);
          if (edgeIter == edgeList->end()) {
            (*edgeList)[pos] = new VariantEdge(pos);
          }

          if (variant1Iter->allele == 0) {
            (*edgeList)[pos]->ref->addSubEdge(variant1Iter->quality,
                                              *variant2Iter,
                                              readIter->read_name, 0, 1);
          } else if (variant1Iter->allele == 1) {
            (*edgeList)[pos]->alt->addSubEdge(variant1Iter->quality,
                                              *variant2Iter,
                                              readIter->read_name, 0, 1);
          }
        }
        ++searchCount;
        ++variant2Iter;
      }
    }
  }
}

void MethylationGraph::connectResults(std::string chrName,
                                      std::vector<int> &passPosition,
                                      bool hasValidSnpData) {
  std::set<int> strongMethylationPoints;
  std::set<int> weakMethylationPoints;
  std::set<int> addedPositions;
  std::vector<int> prepassPosition;

  if (!hasValidSnpData) {
    collectAllModPoints(strongMethylationPoints);
  } else {
    collectStrongPoints(strongMethylationPoints, weakMethylationPoints,
                        prepassPosition);
  }

  expandStrongPoints(strongMethylationPoints, prepassPosition,
                     weakMethylationPoints, addedPositions, hasValidSnpData);
  iterateWeakPoints(weakMethylationPoints, addedPositions, prepassPosition);
  filterByNeighborConnection(prepassPosition, passPosition);
}

float MethylationGraph::computeMinConnection(int pos1, int pos2) const {
  return std::max(((*nodeInfo)[pos1].size() + (*nodeInfo)[pos2].size()) / 4.0,
                  6.0);
}

double MethylationGraph::computeMajorRatio(std::pair<int, int> counts) const {
  return (double)std::max(counts.first, counts.second) /
         (double)(counts.first + counts.second);
}

bool MethylationGraph::isGoodConnection(std::pair<int, int> counts,
                                        float minConn) const {
  int total = counts.first + counts.second;
  return computeMajorRatio(counts) >= params->connectConfidence &&
         total > minConn;
}

void MethylationGraph::collectAllModPoints(std::set<int> &strongPoints) {
  for (auto nodeIter = nodeInfo->begin(); nodeIter != nodeInfo->end();
       ++nodeIter) {
    if (checkVariantType(nodeIter->first) == VariantType::MOD) {
      strongPoints.insert(nodeIter->first);
    }
  }
}

void MethylationGraph::collectStrongPoints(std::set<int> &strongPoints,
                                           std::set<int> &weakPoints,
                                           std::vector<int> &prepass) {
  std::vector<int> hasConnect;

  for (auto nodeIter = nodeInfo->begin(); nodeIter != nodeInfo->end();
       ++nodeIter) {
    int currPos = nodeIter->first;
    auto nextNodeIter = std::next(nodeIter, 1);
    int searchCount = 0;
    if (nextNodeIter == nodeInfo->end()) {
      break;
    }

    auto edgeIter = edgeList->find(currPos);
    if (edgeIter == edgeList->end())
      continue;

    if (checkVariantType(currPos) == VariantType::MOD) {
      auto searchNodeIter = nextNodeIter;
      while (searchNodeIter != nodeInfo->end() &&
             searchCount < params->connectAdjacent) {
        std::pair<int, int> tmp =
            edgeIter->second->findNumberOfRead(searchNodeIter->first);
        int totalConnectReads = tmp.first + tmp.second;
        float minimumConnection =
            computeMinConnection(currPos, searchNodeIter->first);
        if (totalConnectReads <= minimumConnection) {
          break;
        }
        if (checkVariantType(searchNodeIter->first) == VariantType::SNP) {
          hasConnect.push_back(currPos);
          if (isGoodConnection(tmp, minimumConnection) &&
              strongPoints.count(currPos) == 0) {
            strongPoints.insert(currPos);
            break;
          }
        }
        ++searchNodeIter;
        ++searchCount;
      }
      if (std::find(hasConnect.begin(), hasConnect.end(), currPos) ==
          hasConnect.end()) {
        weakPoints.insert(currPos);
      }
    } else if (checkVariantType(currPos) == VariantType::SNP) {
      auto searchNodeIter = nextNodeIter;
      prepass.push_back(currPos);
      while (searchNodeIter != nodeInfo->end()) {
        std::pair<int, int> tmp =
            edgeIter->second->findNumberOfRead(searchNodeIter->first);
        int totalConnectReads = tmp.first + tmp.second;
        float minimumConnection =
            computeMinConnection(currPos, searchNodeIter->first);
        if (totalConnectReads <= minimumConnection) {
          break;
        }
        if (checkVariantType(searchNodeIter->first) == VariantType::MOD) {
          hasConnect.push_back(searchNodeIter->first);
          if (isGoodConnection(tmp, minimumConnection) &&
              strongPoints.count(nextNodeIter->first) == 0) {
            strongPoints.insert(nextNodeIter->first);
          }
        }
        ++searchNodeIter;
        ++searchCount;
      }
    }
  }
}

void MethylationGraph::expandStrongPoints(const std::set<int> &strongPoints,
                                          std::vector<int> &prepass,
                                          std::set<int> &weakPoints,
                                          std::set<int> &addedPositions,
                                          bool hasValidSnpData) {
  for (auto it1 = strongPoints.begin(); it1 != strongPoints.end(); ++it1) {
    int pos1 = *it1;
    auto it2 = std::next(it1, 1);
    auto searchNodeIter = it2;
    int searchCount = 0;
    auto edgeIter = edgeList->find(pos1);
    if (edgeIter == edgeList->end())
      continue;

    while (searchNodeIter != strongPoints.end() &&
           searchCount < params->connectAdjacent) {
      int pos2 = *searchNodeIter;
      std::pair<int, int> tmp = edgeIter->second->findNumberOfRead(pos2);
      int totalConnectReads = tmp.first + tmp.second;
      float minimumConnection = computeMinConnection(pos1, pos2);

      if (totalConnectReads <= minimumConnection) {
        break;
      }

      if (isGoodConnection(tmp, minimumConnection)) {
        if (addedPositions.count(pos1) == 0) {
          prepass.push_back(pos1);
          addedPositions.insert(pos1);
          // Only add positions to weakMethylationPoints if there is valid SNP
          // data (for third pass)
          if (hasValidSnpData) {
            weakPoints.insert(pos1);
          }
        }
        if (addedPositions.count(pos2) == 0) {
          prepass.push_back(pos2);
          addedPositions.insert(pos2);
          // Only add positions to weakMethylationPoints if there is valid SNP
          // data (for third pass)
          if (hasValidSnpData) {
            weakPoints.insert(pos2);
          }
        }
      }
      ++searchNodeIter;
      ++searchCount;
    }
  }
}

void MethylationGraph::iterateWeakPoints(std::set<int> &weakPoints,
                                         std::set<int> &addedPositions,
                                         std::vector<int> &prepass) {
  std::set<int> weakPoints2;
  std::set<int> addedPositions2;

  // third pass: evaluate connections between weak methylation points
  for (int i = 0; i < params->iterCount; i++) {
    // Use alternating sets for each iteration
    auto &currentWeakPoints = (i % 2 == 0) ? weakPoints : weakPoints2;
    auto &nextWeakPoints = (i % 2 == 0) ? weakPoints2 : weakPoints;
    auto &currentAddedPositions =
        (i % 2 == 0) ? addedPositions : addedPositions2;
    auto &nextAddedPositions = (i % 2 == 0) ? addedPositions2 : addedPositions;

    // Clear the target set for this iteration
    nextWeakPoints.clear();
    nextAddedPositions.clear();

    for (auto it1 = currentWeakPoints.begin(); it1 != currentWeakPoints.end();
         ++it1) {
      int currPos = *it1;
      auto nextIter = it1;
      int nextSearchCount = 0;
      bool isAdded = false;
      auto edgeIter = edgeList->find(currPos);
      if (edgeIter == edgeList->end())
        continue;
      while (++nextIter != currentWeakPoints.end() &&
             nextSearchCount < params->connectAdjacent) {
        int nextPos = *nextIter;
        if (currentAddedPositions.count(nextPos) == 0 &&
            currentAddedPositions.count(currPos) == 0) {
          ++nextSearchCount;
          continue;
        }
        isAdded = true;
        std::pair<int, int> tmp = edgeIter->second->findNumberOfRead(nextPos);
        int totalConnectReads = tmp.first + tmp.second;
        float minimumConnection = computeMinConnection(currPos, nextPos);
        if (totalConnectReads <= minimumConnection) {
          break;
        }
        if (isGoodConnection(tmp, minimumConnection)) {
          if (std::find(prepass.begin(), prepass.end(), currPos) ==
              prepass.end()) {
            prepass.push_back(currPos);
            nextWeakPoints.insert(currPos);
            nextAddedPositions.insert(currPos);
          }
          if (std::find(prepass.begin(), prepass.end(), nextPos) ==
              prepass.end()) {
            prepass.push_back(nextPos);
            nextWeakPoints.insert(nextPos);
            nextAddedPositions.insert(nextPos);
          }
        }
        ++nextSearchCount;
      }
      if (!isAdded) {
        nextWeakPoints.insert(currPos);
      }
    }
  }
}

void MethylationGraph::filterByNeighborConnection(
    std::vector<int> &prepass, std::vector<int> &passPosition) {
  // Ensure passPosition is sorted by position
  std::sort(prepass.begin(), prepass.end());
  // Fourth step: Filter positions that do not have good connections to both
  // neighbors
  for (size_t i = 0; i < prepass.size(); ++i) {
    int pos = prepass[i];
    bool hasGoodPrevConnection = false;
    bool hasGoodNextConnection = false;
    if (nodeInfo->find(pos) != nodeInfo->end()) {
      if (checkVariantType(pos) == VariantType::SNP) {
        continue;
      }
    }

    // Check connection with previous position
    if (i > 0) {
      int prevPos = prepass[i - 1];
      auto edgeIter = edgeList->find(prevPos);
      if (edgeIter == edgeList->end()) {
        hasGoodPrevConnection = true;
        continue;
      }
      std::pair<int, int> tmp = edgeIter->second->findNumberOfRead(pos);
      int totalConnectReads = tmp.first + tmp.second;
      double majorRatio = computeMajorRatio(tmp);
      if (totalConnectReads != 0) {

        if (majorRatio >= params->connectConfidence && totalConnectReads >= 6) {
          hasGoodPrevConnection = true;
        }
      }
    }

    // Check connection with next position
    if (i < prepass.size() - 1 && hasGoodPrevConnection) {
      int nextPos = prepass[i + 1];
      auto edgeIter = edgeList->find(pos);
      if (edgeIter == edgeList->end()) {
        hasGoodNextConnection = true;
        continue;
      }
      std::pair<int, int> tmp = edgeIter->second->findNumberOfRead(nextPos);
      int totalConnectReads = tmp.first + tmp.second;
      double majorRatio = computeMajorRatio(tmp);
      if (totalConnectReads != 0) {
        if (majorRatio >= params->connectConfidence && totalConnectReads >= 6) {
          hasGoodNextConnection = true;
        }
      }
    }

    // Only keep positions with good connections to both neighbors (or edge
    // positions)
    if (hasGoodNextConnection || i == 0 || i == prepass.size() - 1) {
      passPosition.push_back(pos);
    }
  }
  prepass.clear();
}

int MethylationGraph::checkVariantType(int position) {
  auto nodeIter = nodeInfo->find(position);
  if (nodeIter != nodeInfo->end()) {
    for (const auto &nodeType : nodeIter->second) {
      if (nodeType.second == VariantType::MOD) {
        return VariantType::MOD; // Return true if any type is MOD
      } else if (nodeType.second == VariantType::SNP) {
        return VariantType::SNP; // Return true if any type is SNP
      } else if (nodeType.second == VariantType::INDEL) {
        return VariantType::INDEL; // Return true if any type is INDEL
      } else if (nodeType.second == VariantType::SV) {
        return VariantType::SV; // Return true if any type is SV
      } else {
        return -1; // Return -1 if no variant type is found
      }
    }
  }
  return -1; // Return -1 if the position is not in the nodeInfo
}

void MethylationGraph::destroy() {

  for (auto edgeIter = edgeList->begin(); edgeIter != edgeList->end();
       edgeIter++) {
    edgeIter->second->ref->destroy();
    edgeIter->second->alt->destroy();
    delete edgeIter->second->ref;
    delete edgeIter->second->alt;
  }

  for (auto nodeIter = forwardModNode->begin();
       nodeIter != forwardModNode->end(); nodeIter++) {
    delete nodeIter->second;
  }

  for (auto nodeIter = reverseModNode->begin();
       nodeIter != reverseModNode->end(); nodeIter++) {
    delete nodeIter->second;
  }

  delete nodeInfo;
  delete edgeList;
}

MethSnpParser::MethSnpParser(ModCallParameters &in_params)
    : hasSnpData(false) {
  chrVariant = new std::map<std::string, std::map<int, SnpVariant>>;
  params = &in_params;

  if (params->snpFile.empty()) {
    std::cerr << "No SNP file provided, running without SNP variants.";
    return;
  }

  hasSnpData = loadSNPVCF();
}

bool MethSnpParser::loadSNPVCF() {
  htsFile *inf = bcf_open(params->snpFile.c_str(), "r");
  if (inf == nullptr) {
    std::cerr << "Warning: Could not open SNP file " << params->snpFile
              << ", running without SNP variants.";
    return false;
  }

  bcf_hdr_t *hdr = bcf_hdr_read(inf);
  if (hdr == nullptr) {
    std::cerr << "Warning: Could not read SNP file header, running without SNP "
                 "variants.";
    bcf_close(inf);
    return false;
  }

  int nseq = 0;
  const char **seqnames = bcf_hdr_seqnames(hdr, &nseq);

  std::string allSmples = "-";
  int is_file = bcf_hdr_set_samples(hdr, allSmples.c_str(), 0);
  if (is_file != 0) {
    std::cout << "error or a positive integer if the list contains samples not "
                 "present in the VCF header\n";
  }

  bcf1_t *rec = bcf_init();
  int ngt_arr = 0;
  int ngt = 0;
  int *gt = NULL;

  while (bcf_read(inf, hdr, rec) == 0) {
    if (bcf_is_snp(rec)) {
      ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);

      if (ngt < 0) {
        std::cerr << "pos " << rec->pos << " missing GT value" << "\n";
        exit(1);
      }

      if ((gt[0] == 2 && gt[1] == 4) || // 0/1
          (gt[0] == 4 && gt[1] == 2) || // 1/0
          (gt[0] == 2 && gt[1] == 5) || // 0|1
          (gt[0] == 4 && gt[1] == 3)    // 1|0
      ) {
        std::string chr = seqnames[rec->rid];
        int variantPos = rec->pos;
        SnpVariant tmp;
        tmp.Ref = rec->d.allele[0];
        tmp.Alt = rec->d.allele[1];

        if (rec->d.allele[1][2] != '\0') {
          continue;
        }

        (*chrVariant)[chr][variantPos] = tmp;
      }
    }
  }

  if (gt)
    free(gt);
  bcf_destroy(rec);
  bcf_hdr_destroy(hdr);
  bcf_close(inf);

  std::cerr << "Successfully loaded SNP variants from " << params->snpFile
            << "\n";
  return true;
}

void MethSnpParser::parserProcess(std::string &input) {}

MethSnpParser::~MethSnpParser() { delete chrVariant; }

bool MethSnpParser::hasValidSnpData() const { return hasSnpData; }

void MethSnpParser::writeDataLine(const std::string &input,
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

  // this pos is phased
  if (psElementIter != phasingResult.end() &&
      posIter != (*chrVariant)[fields[0]].end()) {
    fields[8] = fields[8] + ":PS";
    fields[9] = fields[9] + ":" + std::to_string((*psElementIter).second.block);

    int modify_start = findFieldValueStart(fields[8], fields[9], "GT");
    fields[9][modify_start] = (*psElementIter).second.RAstatus[0];
    fields[9][modify_start + 1] = '|';
    fields[9][modify_start + 2] = (*psElementIter).second.RAstatus[2];
  } else {
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

std::map<int, SnpVariant>
MethSnpParser::getVariants_markindel(std::string chrName,
                                     const std::string &ref) {
  std::map<int, SnpVariant> targetVariants;
  std::map<std::string, std::map<int, SnpVariant>>::iterator chrIter =
      chrVariant->find(chrName);

  if (chrIter != chrVariant->end()) {
    for (auto innerIter = chrIter->second.begin();
         innerIter != chrIter->second.end(); ++innerIter) {
      innerIter->second.is_danger = isDangerIndel(
          innerIter->first, ref, innerIter->second.Ref.length(),
          innerIter->second.Alt.length());
    }
    targetVariants = chrIter->second;
  }

  return targetVariants;
}

void MethSnpParser::writeResult(PhasingResult phasingResult) {
  dispatchWriteResult(params->snpFile, params->resultPrefix + ".vcf",
                      phasingResult);
}
