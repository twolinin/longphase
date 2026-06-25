#include "HaplotagProcess.h"

namespace {

enum class VariantInputKind {
    Compressed,
    Uncompressed,
    Unsupported
};

VariantInputKind classifyVariantInput(const std::string &variantFile)
{
    if( variantFile.find("gz") != std::string::npos ){
        return VariantInputKind::Compressed;
    }
    if( variantFile.find("vcf") != std::string::npos ){
        return VariantInputKind::Uncompressed;
    }
    return VariantInputKind::Unsupported;
}

template <typename Callback>
void readGzLines(const std::string &variantFile, Callback callback)
{
    gzFile file = gzopen(variantFile.c_str(), "rb");
    if(!file){
        std::cerr<< "Fail to open vcf: " << variantFile << "\n";
        return;
    }

    int buffer_size = 1048576; // 1M
    char* buffer = (char*) malloc(buffer_size);
    if(!buffer){
        std::cerr<<"Failed to allocate buffer\n";
        exit(EXIT_FAILURE);
    }
    char* offset = buffer;
        
    while(true) {
        int len = buffer_size - (offset - buffer);
        if (len == 0){
            buffer_size *= 2; // Double the buffer size
            char* new_buffer = (char*) realloc(buffer, buffer_size);
            if(!new_buffer){
                std::cerr<<"Failed to allocate buffer\n";
                free(buffer);
                exit(EXIT_FAILURE);
            }
            buffer = new_buffer;
            offset = buffer + buffer_size / 2; // Update the offset pointer to the end of the old buffer
            len = buffer_size - (offset - buffer);
        }

        len = gzread(file, offset, len);
        if (len == 0) break;    
        if (len <  0){ 
            int err;
            fprintf (stderr, "Error: %s.\n", gzerror(file, &err));
            exit(EXIT_FAILURE);
        }

        char* cur = buffer;
        char* end = offset+len;
        for (char* eol; (cur<end) && (eol = std::find(cur, end, '\n')) < end; cur = eol + 1)
        {
            std::string input = std::string(cur, eol);
            callback(input);
        }
        // any trailing data in [eol, end) now is a partial line
        offset = std::copy(cur, end, buffer);
    }
    gzclose (file);
    free(buffer);
}

template <typename Callback>
void readPlainLines(const std::string &variantFile, Callback callback)
{
    std::ifstream file(variantFile);
    if (!file.is_open()) {
        std::cerr << "Fail to open vcf: " << variantFile << "\n";
        exit(1);
    }
    std::string input;
    while (!file.eof()) {           // preserve original eof() anti-pattern; line emission sequence must remain identical
        std::getline(file, input);
        callback(input);
    }
}

// Extracts the colon-delimited value starting at fieldStart in sample.
std::string extractFieldValue(const std::string &sample, int fieldStart)
{
    size_t next_colon = sample.find(':', (size_t)fieldStart + 1);
    if (next_colon != std::string::npos)
        return sample.substr((size_t)fieldStart, next_colon - (size_t)fieldStart);
    return sample.substr((size_t)fieldStart);
}

// Returns 1 (0|1: ALT on HP2), 0 (1|0: ALT on HP1), or -1 (other).
int resolveHaplotypeDirection(char gt0, char gt2)
{
    if (gt0 == '0' && gt2 == '1') return 1;
    if (gt0 == '1' && gt2 == '0') return 0;
    return -1;
}

std::pair<int, int> scoreIndelSupport(bool readHasIndel, int hp1Length, int hp2Length)
{
    if (readHasIndel) {
        if (hp1Length != 1 && hp2Length == 1)
            return std::make_pair(1, 0);
        if (hp1Length == 1 && hp2Length != 1)
            return std::make_pair(0, 1);
    }
    else {
        if (hp1Length != 1 && hp2Length == 1)
            return std::make_pair(0, 1);
        if (hp1Length == 1 && hp2Length != 1)
            return std::make_pair(1, 0);
    }
    return std::make_pair(0, 0);
}

}

void HaplotagProcess::variantParser(std::string variantFile){

    VariantInputKind inputKind = classifyVariantInput(variantFile);

    switch(inputKind){
        case VariantInputKind::Compressed:
            compressParser(variantFile);
            break;
        case VariantInputKind::Uncompressed:
            unCompressParser(variantFile);
            break;
        case VariantInputKind::Unsupported:
            std::cerr<<"file: "<< variantFile << "\nnot vcf file. please check filename extension\n";
            exit(EXIT_FAILURE);
    }
}

void HaplotagProcess::compressParser(std::string &variantFile){
    if(variantFile=="")
        return;

    readGzLines(variantFile, [this](std::string &input) {
        parserProcess(input);
    });
}

void HaplotagProcess::unCompressParser(std::string &variantFile){
    if (variantFile == "") {
        return;
    }
    readPlainLines(variantFile, [this](std::string &line) {
        parserProcess(line);
    });
}

void HaplotagProcess::parserProcess(std::string &input){
    if (input.substr(0, 2) == "##") {
        parseHeaderLine(input);
        return;
    }
    if (input.substr(0, 1) == "#") {
        return;
    }
    parseRecordLine(input);
}

void HaplotagProcess::parseHeaderLine(const std::string &line){
    if (!parseSnpFile)
        return;
    if (line.find("contig=") != std::string::npos) {
        int id_start  = line.find("ID=")+3;
        int id_end    = line.find(",length=");
        int len_start = id_end+8;
        int len_end   = line.find(">");

        std::string chr = line.substr(id_start, id_end-id_start);
        int chrLen = std::stoi(line.substr(len_start, len_end-len_start));

        chrVec.push_back(chr);
        chrLength[chr] = chrLen;
    }
    if (line.substr(0, 16) == "##FORMAT=<ID=PS,") {
        if (line.find("Type=Integer") != std::string::npos) {
            integerPS = true;
        }
        else if (line.find("Type=String") != std::string::npos) {
            integerPS = false;
            std::cerr << "PS type is String. Auto index to integer ... ";
        }
        else {
            std::cerr << "ERROR: not found PS type (Type=Integer or Type=String).\n";
            exit(EXIT_SUCCESS);
        }
    }
}

void HaplotagProcess::parseRecordLine(const std::string &line){
    std::istringstream iss(line);
    std::vector<std::string> fields((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());

    if (fields.size() == 0)
        return;

    // trans to 0-base
    int pos = std::stoi(fields[1]) - 1;

    // find GT value start (using helper)
    int modifu_start = findFieldValueStart(fields[8], fields[9], "GT");

    // hetero GT
    if ((fields[9][modifu_start] != fields[9][modifu_start+2]) && fields[9][modifu_start+1] == '|') {
        int ps_start = findFieldValueStart(fields[8], fields[9], "PS");
        std::string psValue = extractFieldValue(fields[9], ps_start);
        int haploDir = resolveHaplotypeDirection(fields[9][modifu_start], fields[9][modifu_start+2]);

        if (parseSnpFile)  processSnpRecord(fields, pos, haploDir, psValue);
        if (parseSVFile)   processSvRecord(fields, haploDir);
        if (parseMODFile)  processModRecord(fields, haploDir);
    }
}

void HaplotagProcess::processSnpRecord(const std::vector<std::string> &fields, int pos, int haploDir, const std::string &psValue){
    std::string chr = fields[0];
    SnpVariant tmp;
    tmp.Ref = fields[3];
    tmp.Alt = fields[4];

    chrVariant[chr][pos] = tmp;

    if (integerPS) {
        chrVariantPS[chr][pos] = std::stoi(psValue);
    }
    else {
        std::map<std::string, int>::iterator psIter = psIndex.find(psValue);
        if (psIter == psIndex.end()) {
            psIndex[psValue] = psIndex.size();
        }
        chrVariantPS[chr][pos] = psIndex[psValue];
    }

    if (haploDir == 1) {
        chrVariantHP1[chr][pos] = fields[3];
        chrVariantHP2[chr][pos] = fields[4];
    }
    else if (haploDir == 0) {
        chrVariantHP1[chr][pos] = fields[4];
        chrVariantHP2[chr][pos] = fields[3];
    }
}

void HaplotagProcess::processSvRecord(const std::vector<std::string> &fields, int haploDir){
    std::string chr = fields[0];
    int start = std::stoi(fields[1]);
    std::string info = fields[7];
    size_t svlenPos = info.find("SVLEN=");
    if (svlenPos != std::string::npos) {
        svlenPos += 6;
        size_t semiPos = info.find(';', svlenPos);
        int svlen = std::stoi(info.substr(svlenPos, semiPos - svlenPos));
        if (haploDir == 1) {
            chrRegions[chr].push_back(std::make_tuple(start, svlen, 1));
        }
        else if (haploDir == 0) {
            chrRegions[chr].push_back(std::make_tuple(start, svlen, 0));
        }
    }
}

void HaplotagProcess::processModRecord(const std::vector<std::string> &fields, int haploDir){
    int read_pos = fields[7].find("MR=");
    read_pos = fields[7].find("=", read_pos);
    read_pos++;

    int next_field = fields[7].find(";", read_pos);
    std::string totalRead = fields[7].substr(read_pos, next_field-read_pos);
    std::stringstream totalReadStream(totalRead);

    std::string read;
    while (std::getline(totalReadStream, read, ',')) {
        auto readIter = readSVHapCount.find(read);
        if (readIter == readSVHapCount.end()) {
            readSVHapCount[read][0] = 0;
            readSVHapCount[read][1] = 0;
        }
        readSVHapCount[read][haploDir]++;
    }
}

std::vector<int> HaplotagProcess::collectLastVariantPos() const {
    std::vector<int> last_pos;
    for (auto chr : chrVec) {
        auto mapIter = chrVariantPS.find(chr);
        if (mapIter != chrVariantPS.end()) {
            auto lastVariantIter = mapIter->second.rbegin();
            if (lastVariantIter != mapIter->second.rend()) {
                last_pos.push_back(lastVariantIter->second);
            } else {
                last_pos.push_back(0);
            }
        } else {
            last_pos.push_back(0);
        }
    }
    return last_pos;
}

std::ofstream* HaplotagProcess::initLogFile(const HaplotagParameters &params) {
    if (!params.writeReadLog)
        return NULL;

    std::ofstream *tagResult = new std::ofstream(params.resultPrefix + ".out");
    if (!tagResult->is_open()) {
        std::cerr << "Fail to open write file: " << params.resultPrefix + ".out" << "\n";
        exit(1);
    }
    (*tagResult) << "##snpFile:"             << params.snpFile             << "\n";
    (*tagResult) << "##svFile:"              << params.svFile              << "\n";
    (*tagResult) << "##bamFile:"             << params.bamFile             << "\n";
    (*tagResult) << "##resultPrefix:"        << params.resultPrefix        << "\n";
    (*tagResult) << "##numThreads:"          << params.numThreads          << "\n";
    (*tagResult) << "##region:"              << params.region              << "\n";
    (*tagResult) << "##qualityThreshold:"    << params.qualityThreshold    << "\n";
    (*tagResult) << "##percentageThreshold:" << params.percentageThreshold << "\n";
    (*tagResult) << "#tagSupplementary:"     << params.tagSupplementary    << "\n";
    (*tagResult) << "#svWindowsize:"         << params.svWindow            << "\n";
    (*tagResult) << "#svThreshold:"          << params.svThreshold         << "\n";
    (*tagResult) << "#Read\t"
                 << "Chr\t"
                 << "ReadStart\t"
                 << "Confidnet(%)\t"
                 << "Haplotype\t"
                 << "PhaseSet\t"
                 << "TotalAllele\t"
                 << "HP1Allele\t"
                 << "HP2Allele\t"
                 << "phasingQuality(PQ)\t"
                 << "(Variant,HP)\t"
                 << "(PhaseSet,Variantcount)\n";
    return tagResult;
}

void HaplotagProcess::tagChromosome(const std::string &chr,
                                    const HaplotagParameters &params,
                                    samFile *in, samFile *out,
                                    bam_hdr_t *bamHdr, hts_idx_t *idx,
                                    bam1_t *aln,
                                    const std::string &chr_reference,
                                    std::ofstream *tagResult) {
    std::time_t begin = time(NULL);
    std::cerr<<"chr: " << chr << " ... " ;

    currentChrVariants = chrVariant[chr];
    firstVariantIter = currentChrVariants.begin();
    std::map<int, SnpVariant>::reverse_iterator last = currentChrVariants.rbegin();

    currentchrRegions = chrRegions[chr];
    firstSVIter = currentchrRegions.begin();

    std::string region = !params.region.empty() ? params.region : chr + ":1-" + std::to_string(chrLength[chr]);
    hts_itr_t *iter = sam_itr_querys(idx, bamHdr, region.c_str());
    int result = 0;
    while ((result = sam_itr_multi_next(in, iter, aln)) >= 0) {
        totalAlignment++;
        int flag = aln->core.flag;

        // Intentionally different from phase/modcall filters: uses
        // qualityThreshold and tagSupplementary; duplicate (0x400) is not
        // filtered here; supplementary handling is config-driven.
        if ( aln->core.qual < params.qualityThreshold ){
            totalUnTagCount++;
        }
        else if( (flag & 0x4) != 0 ){
            totalUnmapped++;
            totalUnTagCount++;
        }
        else if( (flag & 0x100) != 0 ){
            totalSecondary++;
            totalUnTagCount++;
        }
        else if( (flag & 0x800) != 0 && params.tagSupplementary == false ){
            totalSupplementary++;
            totalUnTagCount++;
        }
        else if(last == currentChrVariants.rend()){
            totalUnTagCount++;
        }
        else if(int(aln->core.pos) <= (*last).first){

            if( (flag & 0x800) != 0 ){
                totalSupplementary++;
            }

            int pqValue = 0;
            int haplotype = judgeHaplotype(*bamHdr, *aln, chr, params.percentageThreshold, tagResult, pqValue, chr_reference, params.svWindow, params.svThreshold);

            initFlag(aln, "HP");
            initFlag(aln, "PS");
            initFlag(aln, "PQ");

            if (haplotype != 0){
                int psValue = chrVariantPS[chr][(*firstVariantIter).first];
                bam_aux_append(aln, "HP", 'i', sizeof(haplotype), (uint8_t*) &haplotype);
                bam_aux_append(aln, "PS", 'i', sizeof(psValue), (uint8_t*) &psValue);
                bam_aux_append(aln, "PQ", 'i', sizeof(pqValue), (uint8_t*) &pqValue);
                totalTagCount++;
            }
            else{
                totalUnTagCount++;
            }
        }
        else{
            totalUnTagCount++;
        }

        result = sam_write1(out, bamHdr, aln);
    }
    hts_itr_destroy(iter);
    std::cerr<< difftime(time(NULL), begin) << "s\n";
}

void HaplotagProcess::writeUnmappedTail(samFile *in, samFile *out,
                                        bam_hdr_t *bamHdr, hts_idx_t *idx) {
    hts_itr_t *iter_unmap = sam_itr_querys(idx, bamHdr, "*");
    if (iter_unmap == NULL) {
        std::cerr << "WARN: Cannot create iterator for '*' (unmapped tail)\n";
    } else {
        bam1_t *aln_unmap = bam_init1();
        int r = 0;
        while ((r = sam_itr_multi_next(in, iter_unmap, aln_unmap)) >= 0) {
            initFlag(aln_unmap, "HP");
            initFlag(aln_unmap, "PS");
            initFlag(aln_unmap, "PQ");

            totalAlignment++;
            totalUnmapped++;
            totalUnTagCount++;

            int wr = sam_write1(out, bamHdr, aln_unmap);
            if (wr < 0) {
                std::cerr << "ERROR: Failed to write unmapped read: "
                          << bam_get_qname(aln_unmap) << "\n";
            }
        }
        bam_destroy1(aln_unmap);
        hts_itr_destroy(iter_unmap);
    }
}

void HaplotagProcess::resolveRegionScope(const HaplotagParameters &params){
    if (params.region.empty())
        return;

    auto colonPos = params.region.find(":");
    std::string regionChr = (colonPos != std::string::npos)
        ? params.region.substr(0, colonPos)
        : params.region;

    auto chrVecIter = std::find(chrVec.begin(), chrVec.end(), regionChr);
    if (chrVecIter != chrVec.end()) {
        chrVec = std::vector<std::string>{regionChr};
    } else {
        std::cerr << "ERROR: Incorrect chromosome for input region\n" << std::endl;
        exit(1);
    }
}

void HaplotagProcess::tagRead(HaplotagParameters &params){

    resolveRegionScope(params);

    // input file management
    std::string openBamFile = params.bamFile;
    // open bam file
    samFile *in = hts_open(openBamFile.c_str(), "r");
    // load reference file
    hts_set_fai_filename(in, params.fastaFile.c_str() );
    // input reader
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    // header add pg tag
    sam_hdr_add_pg(bamHdr, "longphase", "VN", params.version.c_str(), "CL", params.command.c_str(), NULL);
    // bam file index
    hts_idx_t *idx = NULL;
    // check input bam file
    if (in == NULL) {
        std::cerr<<"ERROR: Cannot open bam file " << openBamFile.c_str() << "\n";
    }
    // check bam file index
    if ((idx = sam_index_load(in, openBamFile.c_str())) == 0) {
        std::cerr<<"ERROR: Cannot open index for bam file " << openBamFile.c_str() << "\n";
    }

    // output file mangement
    std::string writeBamFile = params.resultPrefix + "." + params.outputFormat;
    // open output bam file
    samFile *out = hts_open(writeBamFile.c_str(), (params.outputFormat == "bam" ? "wb" : "wc" ));
    // load reference file
    hts_set_fai_filename(out, params.fastaFile.c_str() );
    // output writer
    int result = sam_hdr_write(out, bamHdr);
    // check index file
    if ((idx = sam_index_load(in, openBamFile.c_str())) == 0) {
        std::cerr<<"ERROR: Cannot open index for bam file\n";
        exit(1);
    }

    std::vector<int> last_pos = collectLastVariantPos();

    // reference fasta parser
    FastaParser fastaParser(params.fastaFile, chrVec, last_pos, params.numThreads);

    std::ofstream *tagResult = initLogFile(params);
    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};
    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(params.numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
    }
    // set thread
    hts_set_opt(in, HTS_OPT_THREAD_POOL, &threadPool);
    hts_set_opt(out, HTS_OPT_THREAD_POOL, &threadPool);
    // initialize an alignment
    bam1_t *aln = bam_init1();

    // loop all chromosome
    for(auto chr : chrVec ){
        tagChromosome(chr, params, in, out, bamHdr, idx, aln,
                      fastaParser.chrString.at(chr), tagResult);
    }

    writeUnmappedTail(in, out, bamHdr, idx);

    if(tagResult!=NULL){
        (*tagResult).close();
    }
    hts_idx_destroy(idx);
    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
    sam_close(in);
    sam_close(out);
    hts_tpool_destroy(threadPool.pool);

    return;
}

void HaplotagProcess::initFlag(bam1_t *aln, std::string flag){

    uint8_t *hpTag = bam_aux_get(aln, flag.c_str() );

    if( hpTag != NULL )
        bam_aux_del(aln, hpTag);

    return;
}

void HaplotagProcess::collectSVWindowEvidence(const bam1_t &aln, int cigarIdx, int n_cigar,
                                              const uint32_t *cigar,
                                              std::vector<std::tuple<int, int, int> >::iterator svIter,
                                              int svWindow, double svThreshold)
{
    if (svIter == currentchrRegions.end())
        return;

    int sv_start = std::get<0>(*svIter);
    int sv_length = std::get<1>(*svIter);
    int sv_end = sv_start + std::abs(sv_length);
    int svHaplotype = std::get<2>(*svIter);
    double sv_region = sv_end - sv_start + 1;

    for (int j = std::max(cigarIdx - svWindow, 0); j < std::min(cigarIdx + svWindow, n_cigar); j++) {
        int cigar_op = bam_cigar_op(cigar[j]);
        int cigar_oplen = bam_cigar_oplen(cigar[j]);
        if (cigar_op == BAM_CDEL && std::abs(sv_region - cigar_oplen) / std::abs(sv_region) < svThreshold) {
            readSVHapCount[bam_get_qname(&aln)][svHaplotype]++;
            break;
        }

        if (cigar_op == BAM_CINS && std::abs(sv_region - cigar_oplen) / std::abs(sv_region) < svThreshold) {
            readSVHapCount[bam_get_qname(&aln)][svHaplotype]++;
            break;
        }
    }
}

void HaplotagProcess::collectMatchOpEvidence(const bam1_t &aln, const std::string &chrName,
                                             int ref_pos, int query_pos, int length,
                                             int cigarIdx, int n_cigar, const uint32_t *cigar,
                                             std::map<int, SnpVariant>::iterator &variantIter,
                                             std::map<int, int> &variantsHP, std::map<int, int> &countPS,
                                             int &hp1Count, int &hp2Count)
{
    while (variantIter != currentChrVariants.end() && (*variantIter).first < ref_pos + length) {

        int refAlleleLen = (*variantIter).second.Ref.length();
        int altAlleleLen = (*variantIter).second.Alt.length();
        int offset = (*variantIter).first - ref_pos;

        if (offset >= 0) {
            uint8_t *q = bam_get_seq(&aln);
            char base_chr = seq_nt16_str[bam_seqi(q, query_pos + offset)];
            std::string base(1, base_chr);

            // currentVariant is SNP
            if (refAlleleLen == 1 && altAlleleLen == 1) {
                // Detected that the base of the read is either REF or ALT.
                if ((base == (*variantIter).second.Ref) || (base == (*variantIter).second.Alt)) {

                    std::map<int, int>::iterator posPSiter = chrVariantPS[chrName].find((*variantIter).first);

                    if (posPSiter == chrVariantPS[chrName].end()) {
                        std::cerr << (*variantIter).first << "\t"
                                  << (*variantIter).second.Ref << "\t"
                                  << (*variantIter).second.Alt << "\n";
                        exit(EXIT_SUCCESS);
                    }
                    else {
                        if (base == chrVariantHP1[chrName][(*variantIter).first]) {
                            hp1Count++;
                            variantsHP[(*variantIter).first]=0;
                        }
                        if (base == chrVariantHP2[chrName][(*variantIter).first]) {
                            hp2Count++;
                            variantsHP[(*variantIter).first]=1;
                        }
                        countPS[chrVariantPS[chrName][(*variantIter).first]]++;
                    }

                }
            }
            // currentVariant is insertion
            else if (refAlleleLen == 1 && altAlleleLen != 1 && cigarIdx + 1 < n_cigar) {

                int hp1Length = chrVariantHP1[chrName][(*variantIter).first].length();
                int hp2Length = chrVariantHP2[chrName][(*variantIter).first].length();
                bool readHasIndel = (ref_pos + length - 1 == (*variantIter).first && bam_cigar_op(cigar[cigarIdx+1]) == 1);
                std::pair<int, int> indelScore = scoreIndelSupport(readHasIndel, hp1Length, hp2Length);

                hp1Count += indelScore.first;
                hp2Count += indelScore.second;
                if (indelScore.first > 0)
                    variantsHP[(*variantIter).first]=0;
                if (indelScore.second > 0)
                    variantsHP[(*variantIter).first]=1;
                countPS[chrVariantPS[chrName][(*variantIter).first]]++;
            }
            // currentVariant is deletion
            else if (refAlleleLen != 1 && altAlleleLen == 1 && cigarIdx + 1 < n_cigar) {

                int hp1Length = chrVariantHP1[chrName][(*variantIter).first].length();
                int hp2Length = chrVariantHP2[chrName][(*variantIter).first].length();
                bool readHasIndel = (ref_pos + length - 1 == (*variantIter).first && bam_cigar_op(cigar[cigarIdx+1]) == 2);
                std::pair<int, int> indelScore = scoreIndelSupport(readHasIndel, hp1Length, hp2Length);

                hp1Count += indelScore.first;
                hp2Count += indelScore.second;
                if (indelScore.first > 0)
                    variantsHP[(*variantIter).first]=0;
                if (indelScore.second > 0)
                    variantsHP[(*variantIter).first]=1;
                countPS[chrVariantPS[chrName][(*variantIter).first]]++;
            }

        }
        variantIter++;
    }
}

void HaplotagProcess::collectDeletionOpEvidence(const bam1_t &aln, const std::string &chrName,
                                                int ref_pos, int query_pos, int length,
                                                const std::string &ref_string,
                                                std::map<int, SnpVariant>::iterator &variantIter,
                                                std::map<int, int> &variantsHP, std::map<int, int> &countPS,
                                                int &hp1Count, int &hp2Count)
{
    if (ref_string == "")
        return;

    int del_len = length;
    if (ref_pos + del_len + 1 == (*variantIter).first) {
        //if( homopolymerLength((*variantIter).first , ref_string) >=3 ){
            // special case
        //}
    }
    else if ((*variantIter).first >= ref_pos && (*variantIter).first < ref_pos + del_len) {
        // check variant in homopolymer
        if (homopolymerLength((*variantIter).first, ref_string) >= 3) {

            int refAlleleLen = (*variantIter).second.Ref.length();
            int altAlleleLen = (*variantIter).second.Alt.length();

            // SNP
            if (refAlleleLen == 1 && altAlleleLen == 1) {
                // get the next match
                char base_chr = seq_nt16_str[bam_seqi(bam_get_seq(&aln), query_pos)];
                std::string base(1, base_chr);

                if (base == chrVariantHP1[chrName][(*variantIter).first]) {
                    hp1Count++;
                    variantsHP[(*variantIter).first]=0;
                }
                if (base == chrVariantHP2[chrName][(*variantIter).first]) {
                    hp2Count++;
                    variantsHP[(*variantIter).first]=1;
                }
                countPS[chrVariantPS[chrName][(*variantIter).first]]++;
            }

            // the read deletion contain VCF's deletion
            else if (refAlleleLen != 1 && altAlleleLen == 1) {

                int hp1Length = chrVariantHP1[chrName][(*variantIter).first].length();
                int hp2Length = chrVariantHP2[chrName][(*variantIter).first].length();
                std::pair<int, int> indelScore = scoreIndelSupport(true, hp1Length, hp2Length);

                hp1Count += indelScore.first;
                hp2Count += indelScore.second;
                if (indelScore.first > 0)
                    variantsHP[(*variantIter).first]=0;
                if (indelScore.second > 0)
                    variantsHP[(*variantIter).first]=1;
                countPS[chrVariantPS[chrName][(*variantIter).first]]++;
            }
        }
    }
}

int HaplotagProcess::decideHaplotype(int hp1Count, int hp2Count, const std::map<int, int> &countPS,
                                     double percentageThreshold, int &pqValue)
{
    double min,max;

    if(hp1Count > hp2Count){
        min = hp2Count;
        max = hp1Count;
    }
    else{
        min = hp1Count;
        max = hp2Count;
    }

    int hpResult = 0;
    if( max/(max+min) < percentageThreshold){
        // no tag
        pqValue = 0;
    }
    else{
        if(hp1Count > hp2Count){
            hpResult = 1;
        }
        if(hp1Count < hp2Count){
            hpResult = 2;
        }
    }

    if( max == 0 ){
        pqValue=0;
    }
    else if( max == ( max + min ) ){
        pqValue=40;
    }
    else{
        pqValue=-10*(std::log10((double)min/double(max+min)));
    }

    // cross two block
    if( countPS.size() > 1  ){
        hpResult = 0;
    }

    return hpResult;
}

void HaplotagProcess::writeReadLogRow(const bam_hdr_t &bamHdr, const bam1_t &aln,
                                      int hpResult, int hp1Count, int hp2Count,
                                      const std::map<int, int> &variantsHP, const std::map<int, int> &countPS,
                                      int pqValue, std::ofstream *tagResult)
{
    if(tagResult==NULL)
        return;

    double min,max;
    if(hp1Count > hp2Count){
        min = hp2Count;
        max = hp1Count;
    }
    else{
        min = hp1Count;
        max = hp2Count;
    }

    //write tag log file
    std::string hpResultStr = ((hpResult == 0 )? "." : std::to_string(hpResult) );
    std::string psResultStr = ".";

    if( hpResultStr != "." ){
        auto psIter = countPS.begin();
        psResultStr = std::to_string((*psIter).first);
    }

    (*tagResult)<< bam_get_qname(&aln)              << "\t"
                << bamHdr.target_name[aln.core.tid] << "\t"
                << aln.core.pos                     << "\t"
                << max/(max+min)                    << "\t"
                << hpResultStr                      << "\t"
                << psResultStr                      << "\t"
                << hp1Count+hp2Count                << "\t"
                << hp1Count                         << "\t"
                << hp2Count                         << "\t"
                << pqValue                          << "\t";


    // print position and HP
    for(auto v : variantsHP ){
        (*tagResult)<< " " << v.first << "," << v.second ;
    }

    (*tagResult) << "\t";

    // belong PS, number of variant
    for(auto v : countPS ){
        (*tagResult)<< " " << v.first << "," << v.second ;
    }

    (*tagResult)<< "\n";
}

int HaplotagProcess::judgeHaplotype(const  bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, const std::string &ref_string, int svWindow, double svThreshold){

    int hp1Count = 0;
    int hp2Count = 0;
    //record variants on this read
    std::map<int,int> variantsHP;
    std::map<int,int> countPS;

    // Skip variants that are to the left of this read
    while( firstVariantIter != currentChrVariants.end() && (*firstVariantIter).first < aln.core.pos ){
        firstVariantIter++;
    }

    // Skip SVs that are to the left of this read
    while (firstSVIter != currentchrRegions.end() && std::get<0>(*firstSVIter) < aln.core.pos)
    {
        firstSVIter++;
    }

    if( firstVariantIter == currentChrVariants.end() ){
        return 0;
    }

    // position relative to reference
    int ref_pos = aln.core.pos;
    // position relative to read
    int query_pos = 0;
    // set variant start for current alignment
    std::map<int, SnpVariant>::iterator currentVariantIter = firstVariantIter;
    // set first SV start for current alignment
    std::vector<std::tuple<int, int, int> >::iterator currentchrRegionIter = firstSVIter;

    // reading cigar to detect snp on this read
    int aln_core_n_cigar = int(aln.core.n_cigar);
    for(int i = 0; i < aln_core_n_cigar ; i++ ){
        uint32_t *cigar = bam_get_cigar(&aln);
        int cigar_op = bam_cigar_op(cigar[i]);
        int length   = bam_cigar_oplen(cigar[i]);

        // iterator next variant
        while( currentVariantIter != currentChrVariants.end() && (*currentVariantIter).first < ref_pos ){
            currentVariantIter++;
        }
         // iterator next SV
        while (currentchrRegionIter != currentchrRegions.end() && std::get<0>(*currentchrRegionIter) < ref_pos)
        {
            currentchrRegionIter++;
        }

        collectSVWindowEvidence(aln, i, aln_core_n_cigar, cigar, currentchrRegionIter, svWindow, svThreshold);

        // CIGAR operators: MIDNSHP=X correspond 012345678
        // 0: alignment match (can be a sequence match or mismatch)
        // 7: sequence match
        // 8: sequence mismatch
        if( cigar_op == 0 || cigar_op == 7 || cigar_op == 8 ){
            collectMatchOpEvidence(aln, chrName, ref_pos, query_pos, length, i, aln_core_n_cigar, cigar,
                                   currentVariantIter, variantsHP, countPS, hp1Count, hp2Count);
            query_pos += length;
            ref_pos += length;
        }
            // 1: insertion to the reference
        else if( cigar_op == 1 ){
            query_pos += length;
        }
            // 2: deletion from the reference
        else if( cigar_op == 2 ){
            collectDeletionOpEvidence(aln, chrName, ref_pos, query_pos, length, ref_string,
                                      currentVariantIter, variantsHP, countPS, hp1Count, hp2Count);
            ref_pos += length;
        }
            // 3: skipped region from the reference
        else if( cigar_op == 3 ){
            ref_pos += length;
        }
            // 4: soft clipping (clipped sequences present in SEQ)
        else if( cigar_op == 4 ){
            query_pos += length;
        }
            // 5: hard clipping (clipped sequences NOT present in SEQ)
            // 6: padding (silent deletion from padded reference)
        else if( cigar_op == 5 || cigar_op == 6 ){
            // do nothing
        }
        else{
            std::cerr<< "alignment find unsupported CIGAR operation from read: " << bam_get_qname(&aln) << "\n";
            exit(1);
        }
    }

    auto readIter = readSVHapCount.find(bam_get_qname(&aln));
    if( readIter != readSVHapCount.end() ){
        hp1Count += readSVHapCount[bam_get_qname(&aln)][0];
        hp2Count += readSVHapCount[bam_get_qname(&aln)][1];
    }

    int hpResult = decideHaplotype(hp1Count, hp2Count, countPS, percentageThreshold, pqValue);

    writeReadLogRow(bamHdr, aln, hpResult, hp1Count, hp2Count, variantsHP, countPS, pqValue, tagResult);

    return hpResult;
}

HaplotagProcess::HaplotagProcess(HaplotagParameters params):
totalAlignment(0),totalSupplementary(0),totalSecondary(0),totalUnmapped(0),totalTagCount(0),totalUnTagCount(0),processBegin(time(NULL)),integerPS(false),parseSnpFile(false),parseSVFile(false),parseMODFile(false)
{
    std::cerr<< "phased SNP file    : " << params.snpFile             << "\n";
    std::cerr<< "phased SV file     : " << params.svFile              << "\n";
    std::cerr<< "phased MOD file    : " << params.modFile             << "\n";
    std::cerr<< "input bam file     : " << params.bamFile             << "\n";
    std::cerr<< "input ref file     : " << params.fastaFile           << "\n";
    std::cerr<< "output bam file    : " << params.resultPrefix + "." + params.outputFormat << "\n";
    std::cerr<< "number of threads  : " << params.numThreads          << "\n";
    std::cerr<< "write log file     : " << (params.writeReadLog ? "true" : "false") << "\n";
    std::cerr<< "log file           : " << (params.writeReadLog ? (params.resultPrefix+".out") : "") << "\n";
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "tag region                   : " << (!params.region.empty() ? params.region : "all") << "\n";
    std::cerr<< "filter mapping quality below : " << params.qualityThreshold    << "\n";
    std::cerr<< "percentage threshold         : " << params.percentageThreshold << "\n";
    if (!params.svFile.empty()) {
        std::cerr<< "SV windowsize                : " << params.svWindow    << "\n";
        std::cerr<< "SV threshold                 : " << params.svThreshold << "\n";
    }
    std::cerr<< "tag supplementary            : " << (params.tagSupplementary ? "true" : "false") << "\n";
    std::cerr<< "-------------------------------------------\n";

    // load SNP vcf file
    std::time_t begin = time(NULL);
    std::cerr<< "parsing SNP VCF ... ";
    parseSnpFile = true;
    variantParser(params.snpFile);
    parseSnpFile = false;
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // load SV vcf file
    if(params.svFile!=""){
        begin = time(NULL);
        std::cerr<< "parsing SV VCF ... ";
        parseSVFile = true;
        variantParser(params.svFile);
        parseSVFile = false;
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }

    // load MOD vcf file
    if(params.modFile!=""){
        begin = time(NULL);
        std::cerr<< "parsing MOD VCF ... ";
        parseMODFile = true;
        variantParser(params.modFile);
        parseMODFile = false;
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }

    // tag read
    begin = time(NULL);
    std::cerr<< "tag read start ...\n";
    tagRead(params);
    std::cerr<< "tag read " << difftime(time(NULL), begin) << "s\n";

    return;
};

HaplotagProcess::~HaplotagProcess(){
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "total process time:  " << difftime(time(NULL), processBegin) << "s\n";
    std::cerr<< "total alignment:     " << totalAlignment     << "\n";
    std::cerr<< "total supplementary: " << totalSupplementary << "\n";
    std::cerr<< "total secondary:     " << totalSecondary     << "\n";
    std::cerr<< "total unmapped:      " << totalUnmapped      << "\n";
    std::cerr<< "total tag alignment: " << totalTagCount     << "\n";
    std::cerr<< "total untagged:      " << totalUnTagCount   << "\n";
};
