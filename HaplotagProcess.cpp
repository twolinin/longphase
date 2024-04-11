#include "HaplotagProcess.h"


void HaplotagProcess::variantParser(std::string variantFile){

    if( variantFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressParser(variantFile);
    }
    else if( variantFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressParser(variantFile);
    }
    else{
        std::cerr<<"file: "<< variantFile << "\nnot vcf file. please check filename extension\n";
        exit(EXIT_FAILURE);
    }
    return;
}

void HaplotagProcess::compressParser(std::string &variantFile){
    gzFile file = gzopen(variantFile.c_str(), "rb");
    if(variantFile=="")
        return;
    if(!file){
        std::cerr<< "Fail to open vcf: " << variantFile << "\n";
    }
    else{  
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
                parserProcess(input);
            }
            // any trailing data in [eol, end) now is a partial line
            offset = std::copy(cur, end, buffer);
        }
        gzclose (file);
        free(buffer);
    }    
}

void HaplotagProcess::unCompressParser(std::string &variantFile){
    std::ifstream originVcf(variantFile);
    if(variantFile=="")
        return;
    if(!originVcf.is_open()){
        std::cerr<< "Fail to open vcf: " << variantFile << "\n";
        exit(1);
    }
    else{
        std::string input;
        while(! originVcf.eof() ){
            std::getline(originVcf, input);
            parserProcess(input);
        }
    }
}

void HaplotagProcess::parserProcess(std::string &input){
    if( input.substr(0, 2) == "##" && parseSnpFile){
        if( input.find("contig=")!= std::string::npos ){
            int id_start  = input.find("ID=")+3;
            int id_end    = input.find(",length=");
            int len_start = id_end+8;
            int len_end   = input.find(">");
            
            std::string chr = input.substr(id_start,id_end-id_start);
            int chrLen = std::stoi( input.substr(len_start,len_end-len_start) );

            chrVec.push_back(chr);
            chrLength[chr]=chrLen;
        }
        if( input.substr(0, 16) == "##FORMAT=<ID=PS," ){
            if( input.find("Type=Integer")!= std::string::npos ){
                integerPS = true;
            }
            else if( input.find("Type=String")!= std::string::npos ){
                integerPS = false;
                std::cerr<< "PS type is String. Auto index to integer ... ";
            }
            else{
                std::cerr<< "ERROR: not found PS type (Type=Integer or Type=String).\n"; 
                exit(EXIT_SUCCESS);
            }
        }
    }
    else if ( input.substr(0, 1) == "#" ){
        
    }
    else{
        std::istringstream iss(input);
        std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

        if( fields.size() == 0 )
            return;
            
        // trans to 0-base
        int pos = std::stoi( fields[1] ) - 1;
        std::string chr = fields[0];
                
        // find GT flag colon
        int colon_pos = 0;
        int gt_pos = fields[8].find("GT");
        for(int i =0 ; i< gt_pos ; i++){
            if(fields[8][i]==':')
                colon_pos++;
        }
        // find GT value start
        int current_colon = 0;
        int modifu_start = 0;
        for(unsigned int i =0; i < fields[9].length() ; i++){
            if( current_colon >= colon_pos )
                break;
            if(fields[9][i]==':')
                current_colon++;  
            modifu_start++;
        }

        // hetero GT
        if( (fields[9][modifu_start] != fields[9][modifu_start+2]) && fields[9][modifu_start+1] == '|' ){
            // find PS flag
            colon_pos = 0;
            int ps_pos = fields[8].find("PS");
            for(int i =0 ; i< ps_pos ; i++){
                if(fields[8][i]==':')
                    colon_pos++;
            }

            // find PS value start
            current_colon = 0;
            int ps_start = 0;
            for(unsigned int i =0; i < fields[9].length() ; i++){
                if( current_colon >= colon_pos )
                    break;
                if(fields[9][i]==':')
                    current_colon++;  
                ps_start++;
            }
            
            std::string psValue;
            // get PS value
            if( fields[9].find(":",ps_start+1) != std::string::npos ){
                int ps_end_pos = fields[9].find(":",ps_start+1);
                psValue = fields[9].substr(ps_start, ps_end_pos - ps_start);
            }
            else{
                psValue = fields[9].substr(ps_start, fields[9].length() - ps_start);
            }
            
            // snp file
            if( parseSnpFile ){
                RefAlt tmp;
                tmp.Ref = fields[3]; 
                tmp.Alt = fields[4];

                chrVariant[chr][pos]=tmp;
                
                if(integerPS){
                    chrVariantPS[chr][pos]=std::stoi(psValue);
                }
                else{
                    std::map<std::string, int>::iterator psIter = psIndex.find(psValue);
                    
                    if( psIter == psIndex.end() ){
                        psIndex[psValue] = psIndex.size();
                    }
                    chrVariantPS[chr][pos]=psIndex[psValue];
                }
                
                // record haplotype allele
                if( fields[9][modifu_start] == '0' && fields[9][modifu_start+2] == '1' ){
                    chrVariantHP1[chr][pos]=fields[3];
                    chrVariantHP2[chr][pos]=fields[4];
                }
                else if( fields[9][modifu_start] == '1' && fields[9][modifu_start+2] == '0' ){
                    chrVariantHP1[chr][pos]=fields[4];
                    chrVariantHP2[chr][pos]=fields[3];
                }
            }
            // sv file
            if( parseSVFile ){
                // get read INFO
                int read_pos = fields[7].find("RNAMES=");
                read_pos = fields[7].find("=",read_pos);
                read_pos++;
                        
                int next_field = fields[7].find(";",read_pos);
                std::string totalRead = fields[7].substr(read_pos,next_field-read_pos);
                std::stringstream totalReadStream(totalRead);
                
                int svHaplotype;
                // In which haplotype does SV occur
                if( fields[9][modifu_start] == '0' && fields[9][modifu_start+2] == '1' ){
                    svHaplotype = 1;
                }
                else if( fields[9][modifu_start] == '1' && fields[9][modifu_start+2] == '0' ){
                    svHaplotype = 0;
                }
                
                std::string read;
                while(std::getline(totalReadStream, read, ','))
                {
                   auto readIter = readSVHapCount.find(read);
                   if(readIter==readSVHapCount.end()){
                       readSVHapCount[read][0]=0;
                       readSVHapCount[read][1]=0;
                   }
                   readSVHapCount[read][svHaplotype]++;
                }
                
            }
            // mod file
            if( parseMODFile ){
                // get read INFO
                int read_pos = fields[7].find("MR=");
                read_pos = fields[7].find("=",read_pos);
                read_pos++;
                        
                int next_field = fields[7].find(";",read_pos);
                std::string totalRead = fields[7].substr(read_pos,next_field-read_pos);
                std::stringstream totalReadStream(totalRead);
                
                int modHaplotype;
                // In which haplotype does SV occur
                if( fields[9][modifu_start] == '0' && fields[9][modifu_start+2] == '1' ){
                    modHaplotype = 1;
                }
                else if( fields[9][modifu_start] == '1' && fields[9][modifu_start+2] == '0' ){
                    modHaplotype = 0;
                }
                
                std::string read;
                while(std::getline(totalReadStream, read, ','))
                {
                   auto readIter = readSVHapCount.find(read);
                   if(readIter==readSVHapCount.end()){
                       readSVHapCount[read][0]=0;
                       readSVHapCount[read][1]=0;
                   }
                   readSVHapCount[read][modHaplotype]++;
                }
            }
        }
    }
}

void HaplotagProcess::tagRead(HaplotagParameters &params){

    // update chromosome processing based on region
    if (!params.region.empty()) {
        auto colonPos = params.region.find(":");
        std::string regionChr;
        if (colonPos != std::string::npos) {
            regionChr = params.region.substr(0, colonPos);
        }
        else {
            regionChr = params.region;
        }
        auto chrVecIter = std::find(chrVec.begin(), chrVec.end(), regionChr);
        if (chrVecIter != chrVec.end()) {
            chrVec = std::vector<std::string>{regionChr};
        } else {
            std::cerr << "ERROR: Incorrect chromosome for input region\n" << std::endl;
            exit(1);
        }
    }
    // input file management
    std::string openBamFile = params.bamFile;
    // open bam file
    samFile *in = hts_open(openBamFile.c_str(), "r");
    // load reference file
    hts_set_fai_filename(in, params.fastaFile.c_str() );
    // input reader
    bam_hdr_t *bamHdr = sam_hdr_read(in);
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
    std::string writeBamFile = params.resultPrefix + ".bam";
    // open output bam file
    samFile *out = hts_open(writeBamFile.c_str(), "wb");
    // output writer
    int result = sam_hdr_write(out, bamHdr);
    // check index file
    if ((idx = sam_index_load(in, openBamFile.c_str())) == 0) {
        std::cerr<<"ERROR: Cannot open index for bam file\n";
        exit(1);
    }

    // record reference last variant pos
    std::vector<int> last_pos;
    for( auto chr : chrVec ){
        auto lastVariantIter = chrVariantPS[chr].rbegin();
        if( lastVariantIter != chrVariantPS[chr].rend() ){
            last_pos.push_back(lastVariantIter->second);
        }
        else{
            last_pos.push_back(0);
        }
    }

    // reference fasta parser
    FastaParser fastaParser(params.fastaFile, chrVec, last_pos, params.numThreads);

    // write tag read detail information
    std::ofstream *tagResult=NULL;
    if(params.writeReadLog){
        tagResult=new std::ofstream(params.resultPrefix+".out");
        if(!tagResult->is_open()){
            std::cerr<< "Fail to open write file: " << params.resultPrefix+".out" << "\n";
            exit(1);
        }
        else{
            (*tagResult) << "##snpFile:"             << params.snpFile             << "\n";
            (*tagResult) << "##svFile:"              << params.svFile              << "\n";
            (*tagResult) << "##bamFile:"             << params.bamFile             << "\n";
            (*tagResult) << "##resultPrefix:"        << params.resultPrefix        << "\n";
            (*tagResult) << "##numThreads:"          << params.numThreads          << "\n";
            (*tagResult) << "##region:"              << params.region              << "\n";
            (*tagResult) << "##qualityThreshold:"    << params.qualityThreshold    << "\n";
            (*tagResult) << "##percentageThreshold:" << params.percentageThreshold << "\n";
            (*tagResult) << "#tagSupplementary:"     << params.tagSupplementary    << "\n";
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

        }
    }
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
        std::time_t begin = time(NULL);
        std::cerr<<"chr: " << chr << " ... " ;
        // records all variants within this chromosome.
        currentChrVariants = chrVariant[chr];
        // since each read is sorted based on the start coordinates, to save time, 
        // firstVariantIter keeps track of the first variant that each read needs to check.
        firstVariantIter = currentChrVariants.begin();
        // get the coordinates of the last variant
        // the tagging process will not be perform if the read's start coordinate are over than last variant.
        std::map<int, RefAlt>::reverse_iterator last = currentChrVariants.rbegin();
        
        // fetch chromosome string
        std::string chr_reference = fastaParser.chrString.at(chr);
        
        // tagging will be attempted for reads within the specified coordinate region.
        std::string region = !params.region.empty() ? params.region : chr + ":1-" + std::to_string(chrLength[chr]);
        hts_itr_t *iter = sam_itr_querys(idx, bamHdr, region.c_str());
        // iter all reads
        while ((result = sam_itr_multi_next(in, iter, aln)) >= 0) {
            totalAlignment++;
            int flag = aln->core.flag;

            if ( aln->core.qual < params.qualityThreshold ){
                // mapping quality is lower than threshold
		totalUnTagCount++;
            }
            else if( (flag & 0x4) != 0 ){
                // read unmapped
                totalUnmapped++;
		totalUnTagCount++;
            }
            else if( (flag & 0x100) != 0 ){
                // secondary alignment. repeat.
                // A secondary alignment occurs when a given read could align reasonably well to more than one place.
                totalSecondary++;
		totalUnTagCount++;
            }
            else if( (flag & 0x800) != 0 && params.tagSupplementary == false ){
                // supplementary alignment
                // A chimeric alignment is represented as a set of linear alignments that do not have large overlaps.
                totalSupplementary++;
		totalUnTagCount++;
            }
            else if(last == currentChrVariants.rend()){
                // skip 
                totalUnTagCount++;
            }
            else if(int(aln->core.pos) <= (*last).first){
                
		if( (flag & 0x800) != 0 ){
		    totalSupplementary++;
		}

		int pqValue = 0;
                int haplotype = judgeHaplotype(*bamHdr, *aln, chr, params.percentageThreshold, tagResult, pqValue, chr_reference);

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
                else {
                    totalUnTagCount++;
                }
            }
	    else {
                totalUnTagCount++;
	    }

            // write this alignment to result bam file
            result = sam_write1(out, bamHdr, aln);
        }
        std::cerr<< difftime(time(NULL), begin) << "s\n";
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

int HaplotagProcess::judgeHaplotype(const  bam_hdr_t &bamHdr,const bam1_t &aln, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue, const std::string &ref_string){

    int hp1Count = 0;
    int hp2Count = 0;
    //record variants on this read
    std::map<int,int> variantsHP;
    std::map<int,int> countPS;

    // Skip variants that are to the left of this read
    while( firstVariantIter != currentChrVariants.end() && (*firstVariantIter).first < aln.core.pos ){
        firstVariantIter++;
    }
    
    if( firstVariantIter == currentChrVariants.end() ){
        return 0;
    }

    // position relative to reference
    int ref_pos = aln.core.pos;
    // position relative to read
    int query_pos = 0;
    // set variant start for current alignment
    std::map<int, RefAlt>::iterator currentVariantIter = firstVariantIter;

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

        // CIGAR operators: MIDNSHP=X correspond 012345678
        // 0: alignment match (can be a sequence match or mismatch)
        // 7: sequence match
        // 8: sequence mismatch
        if( cigar_op == 0 || cigar_op == 7 || cigar_op == 8 ){

            while( currentVariantIter != currentChrVariants.end() && (*currentVariantIter).first < ref_pos + length){

                int refAlleleLen = (*currentVariantIter).second.Ref.length();
                int altAlleleLen = (*currentVariantIter).second.Alt.length();
                int offset = (*currentVariantIter).first - ref_pos;

                if( offset < 0){
                }
                else{
                    uint8_t *q = bam_get_seq(&aln);
                    char base_chr = seq_nt16_str[bam_seqi(q,query_pos + offset)];
                    std::string base(1, base_chr);
                    
                    // currentVariant is SNP
                    if( refAlleleLen == 1 && altAlleleLen == 1 ){
                        // Detected that the base of the read is either REF or ALT. 
                        if( (base == (*currentVariantIter).second.Ref) || (base == (*currentVariantIter).second.Alt) ){

                            std::map<int, int>::iterator posPSiter = chrVariantPS[chrName].find((*currentVariantIter).first);

                            if( posPSiter == chrVariantPS[chrName].end() ){
                                std::cerr<< (*currentVariantIter).first << "\t"
                                         << (*currentVariantIter).second.Ref << "\t"
                                         << (*currentVariantIter).second.Alt << "\n";
                                exit(EXIT_SUCCESS);
                            }
                            else{
                                if( base == chrVariantHP1[chrName][(*currentVariantIter).first]){
                                    hp1Count++;
                                    variantsHP[(*currentVariantIter).first]=0;
                                }
                                if( base == chrVariantHP2[chrName][(*currentVariantIter).first]){
                                    hp2Count++;
                                    variantsHP[(*currentVariantIter).first]=1;
                                }
                                countPS[chrVariantPS[chrName][(*currentVariantIter).first]]++;
                            }
                            
                        }
                    }
                    // currentVariant is insertion
                    else if( refAlleleLen == 1 && altAlleleLen != 1 && i+1 < aln_core_n_cigar){
                        
                        int hp1Length = chrVariantHP1[chrName][(*currentVariantIter).first].length();
                        int hp2Length = chrVariantHP2[chrName][(*currentVariantIter).first].length();
                        
                        if ( ref_pos + length - 1 == (*currentVariantIter).first && bam_cigar_op(cigar[i+1]) == 1 ) {
                            // hp1 occur insertion
                            if( hp1Length != 1 && hp2Length == 1 ){
                                hp1Count++;
                                variantsHP[(*currentVariantIter).first]=0;
                            }
                            // hp2 occur insertion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp2Count++;
                                variantsHP[(*currentVariantIter).first]=1;
                            }
                        }
                        else {
                            // hp1 occur insertion
                            if( hp1Length != 1 && hp2Length == 1 ){
                                hp2Count++;
                                variantsHP[(*currentVariantIter).first]=1;
                            }
                            // hp2 occur insertion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp1Count++;
                                variantsHP[(*currentVariantIter).first]=0;
                            }
                        }
                        countPS[chrVariantPS[chrName][(*currentVariantIter).first]]++;
                    } 
                    // currentVariant is deletion
                    else if( refAlleleLen != 1 && altAlleleLen == 1 && i+1 < aln_core_n_cigar) {
                        
                        int hp1Length = chrVariantHP1[chrName][(*currentVariantIter).first].length();
                        int hp2Length = chrVariantHP2[chrName][(*currentVariantIter).first].length();
                        
                        if ( ref_pos + length - 1 == (*currentVariantIter).first && bam_cigar_op(cigar[i+1]) == 2 ) {
                            // hp1 occur deletion
                            if( hp1Length != 1 && hp2Length == 1 ){
                                hp1Count++;
                                variantsHP[(*currentVariantIter).first]=0;
                            }
                            // hp2 occur deletion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp2Count++;
                                variantsHP[(*currentVariantIter).first]=1;
                            }
                        }
                        else {
                            // hp2 occur deletion
                            if( hp1Length != 1 && hp2Length == 1 ){
                                hp2Count++;
                                variantsHP[(*currentVariantIter).first]=1;
                            }
                            // hp1 occur deletion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp1Count++;
                                variantsHP[(*currentVariantIter).first]=0;
                            }
                        }
                        countPS[chrVariantPS[chrName][(*currentVariantIter).first]]++;
                    } 

                }
                currentVariantIter++;
            }
            query_pos += length;
            ref_pos += length;
        }
            // 1: insertion to the reference
        else if( cigar_op == 1 ){
            query_pos += length;
        }
            // 2: deletion from the reference
        else if( cigar_op == 2 ){
            
            if(ref_string != ""){
                int del_len = length;
                if ( ref_pos + del_len + 1 == (*currentVariantIter).first ){
                    //if( homopolymerLength((*currentVariantIter).first , ref_string) >=3 ){
                        // special case
                    //}
                }
                else if( (*currentVariantIter).first >= ref_pos  && (*currentVariantIter).first < ref_pos + del_len ){
                    // check variant in homopolymer
                    if( homopolymerLength((*currentVariantIter).first , ref_string) >=3 ){
                        
                        int refAlleleLen = (*currentVariantIter).second.Ref.length();
                        int altAlleleLen = (*currentVariantIter).second.Alt.length();
                        
                        // SNP
                        if( refAlleleLen == 1 && altAlleleLen == 1){
                            // get the next match
                            char base_chr = seq_nt16_str[bam_seqi(bam_get_seq(&aln), query_pos)];
                            std::string base(1, base_chr);

                            if( base == chrVariantHP1[chrName][(*currentVariantIter).first]){
                                hp1Count++;
                                variantsHP[(*currentVariantIter).first]=0;
                            }
                            if( base == chrVariantHP2[chrName][(*currentVariantIter).first]){
                                hp2Count++;
                                variantsHP[(*currentVariantIter).first]=1;
                            }
                            countPS[chrVariantPS[chrName][(*currentVariantIter).first]]++;
                        }
                        
                        // the read deletion contain VCF's deletion
                        else if( refAlleleLen != 1 && altAlleleLen == 1){

                            int hp1Length = chrVariantHP1[chrName][(*currentVariantIter).first].length();
                            int hp2Length = chrVariantHP2[chrName][(*currentVariantIter).first].length();
                            // hp1 occur deletion
                            if( hp1Length != 1 && hp2Length == 1 ){
                                hp1Count++;
                                variantsHP[(*currentVariantIter).first]=0;
                            }
                            // hp2 occur deletion
                            else if( hp1Length == 1 && hp2Length != 1 ){
                                hp2Count++;
                                variantsHP[(*currentVariantIter).first]=1;
                            }
                            countPS[chrVariantPS[chrName][(*currentVariantIter).first]]++;
                        }
                    }
                }
            }
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

    double min,max;

    auto readIter = readSVHapCount.find(bam_get_qname(&aln));
    if( readIter != readSVHapCount.end() ){
        hp1Count += readSVHapCount[bam_get_qname(&aln)][0];
        hp2Count += readSVHapCount[bam_get_qname(&aln)][1];
    }

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

    if(tagResult!=NULL){
        //write tag log file
        std::string hpResultStr = ((hpResult == 0 )? "." : std::to_string(hpResult) );
        std::string psResultStr = ".";

        if( hpResultStr != "." ){
            auto psIter = countPS.begin();
            psResultStr = std::to_string((*psIter).first);
        }

        // cross two block
        if( countPS.size() > 1  ){
            hpResultStr = ".";
            psResultStr = ".";
            hpResult = 0;
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

    return hpResult;
}

HaplotagProcess::HaplotagProcess(HaplotagParameters params):
totalAlignment(0),totalSupplementary(0),totalSecondary(0),totalUnmapped(0),totalTagCount(0),totalUnTagCount(0),processBegin(time(NULL)),integerPS(false)
{
    std::cerr<< "phased SNP file:   " << params.snpFile             << "\n";
    std::cerr<< "phased SV file:    " << params.svFile              << "\n";
    std::cerr<< "phased MOD file:   " << params.modFile              << "\n";
    std::cerr<< "input bam file:    " << params.bamFile             << "\n";
    std::cerr<< "input ref file:    " << params.fastaFile           << "\n";
    std::cerr<< "output bam file:   " << params.resultPrefix+".bam" << "\n";
    std::cerr<< "number of threads: " << params.numThreads          << "\n";
    std::cerr<< "write log file:    " << (params.writeReadLog ? "true" : "false") << "\n";
    std::cerr<< "log file:          " << (params.writeReadLog ? (params.resultPrefix+".out") : "") << "\n";
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "tag region:                    " << (!params.region.empty() ? params.region : "all") << "\n";
    std::cerr<< "filter mapping quality below:  " << params.qualityThreshold    << "\n";
    std::cerr<< "percentage threshold:          " << params.percentageThreshold << "\n";
    std::cerr<< "tag supplementary:             " << (params.tagSupplementary ? "true" : "false") << "\n";
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

