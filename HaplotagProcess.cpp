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
        char buffer[1048576]; // 1M
        char* offset = buffer;
            
        while(true) {
            int len = sizeof(buffer)-(offset-buffer);
            if (len == 0){
                std::cerr<<"Buffer to small for input line lengths\n";
                exit(EXIT_FAILURE);
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
    if( input.find("#")!= std::string::npos && !parseSVFile ){
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
        if( input.find("ID=PS")!= std::string::npos ){
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
    else if( input.find("#")== std::string::npos ){
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
            if( !parseSVFile ){
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
        }
    }
}

void HaplotagProcess::tagRead(HaplotagParameters &params){
    
    // input file management
    std::string openBamFile = params.bamFile;
    // open bam file
    samFile *in = hts_open(openBamFile.c_str(), "r");
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
            (*tagResult) << "##qualityThreshold:"    << params.qualityThreshold    << "\n";
            (*tagResult) << "##percentageThreshold:" << params.percentageThreshold << "\n";
            (*tagResult) << "#tagSupplementary:"     << params.tagSupplementary    << "\n";
            (*tagResult) << "#Read\t" 
                      << "Chr\t"
                      << "ReadStart\t"
                      << "Confidnet(%)\t"
                      << "Haplotype\t"
                      << "TotalAllele\t"
                      << "HP1Allele\t"
                      << "HP2Allele\t"
                      << "phasingQuality(PQ)\t"
                      << "(Variant,HP)\n";
                      
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

        int chunkSize = chrLength[chr]/params.numThreads + params.numThreads;
        char *bamList[params.numThreads];

        currentVariants = chrVariant[chr];
        firstVariantIter = currentVariants.begin();

        std::map<int, RefAlt>::reverse_iterator last = currentVariants.rbegin();

        for(int i=0;i<params.numThreads;i++){
            std::string tmp = chr + ":" + std::to_string(i*chunkSize+1) + "-" + std::to_string((i+1)*chunkSize);
            bamList[i] = new char[tmp.length()+1];
            strcpy(bamList[i], tmp.c_str());
        }
        
        hts_itr_multi_t *iter;
        if( (iter = sam_itr_regarray(idx, bamHdr, bamList, params.numThreads)) == 0){
            std::cerr<<"Warning: Cannot open iterator for " << chr << " for bam file\n";
            continue;
        }  

        while ((result = sam_itr_multi_next(in, iter, aln)) >= 0) {
            totalAlignment++;
            int flag = aln->core.flag;
            
            if ( aln->core.qual < params.qualityThreshold ){
                // mapping quality is lower than threshold
                result = sam_write1(out, bamHdr, aln);
                continue;
            }
            
            if( (flag & 0x4) != 0 ){  
                // read unmapped
                totalUnmapped++;
                result = sam_write1(out, bamHdr, aln);
                continue;
            }
            if( (flag & 0x100) != 0 ){
                // secondary alignment. repeat. 
                // A secondary alignment occurs when a given read could align reasonably well to more than one place.
                totalSecondary++;
                result = sam_write1(out, bamHdr, aln);
                continue;
            }
            if( (flag & 0x800) != 0 && params.tagSupplementary == false ){  
                // supplementary alignment
                // A chimeric alignment is represented as a set of linear alignments that do not have large overlaps.
                totalSupplementary++;
                result = sam_write1(out, bamHdr, aln);
                continue;
            }
            
            Alignment tmp;
            
            tmp.chr = bamHdr->target_name[aln->core.tid];
            tmp.refStart = aln->core.pos;
            tmp.qname = bam_get_qname(aln);
            tmp.qlen = aln->core.l_qseq;
            tmp.cigar_len = aln->core.n_cigar;
            tmp.op = (int*)malloc(tmp.cigar_len*sizeof(int));
            tmp.ol = (int*)malloc(tmp.cigar_len*sizeof(int));
            tmp.is_reverse = bam_is_rev(aln);
            
            //quality string
            uint8_t *q = bam_get_seq(aln); 
            //length of the read
            uint32_t len = aln->core.l_qseq; 
            // set string size
            tmp.qseq = (char *)malloc(len+1);
            // gets nucleotide id and converts them into IUPAC id.
            for(unsigned int i=0; i< len ; i++){
                tmp.qseq[i] = seq_nt16_str[bam_seqi(q,i)]; 
            }

            // set string size
            tmp.quality = (char *)malloc(len+1);
            uint8_t *quality = bam_get_qual(aln);
            // get base quality
            for(unsigned int i=0; i< len ; i++){
                tmp.quality[i] = quality[i];
            }
     
            uint32_t *cigar = bam_get_cigar(aln);
            // store cigar
            for(unsigned int k =0 ; k < aln->core.n_cigar ; k++){
                tmp.op[k] = bam_cigar_op(cigar[k]);
                tmp.ol[k] = bam_cigar_oplen(cigar[k]);
            }
            
            if(last == currentVariants.rend()){
                totalUnTagCuonnt++;
            }
            else if(tmp.refStart <= (*last).first){
                int pqValue = 0;
                int haplotype = judgeHaplotype(tmp, chr, params.percentageThreshold, tagResult, pqValue);
                
                initFlag(aln, "HP");
                initFlag(aln, "PS");
                initFlag(aln, "PQ");
                
                if (haplotype != 0) {
                    
                    int psValue = chrVariantPS[chr][(*firstVariantIter).first];
                    bam_aux_append(aln, "HP", 'i', sizeof(haplotype), (uint8_t*) &haplotype);
                    bam_aux_append(aln, "PS", 'i', sizeof(psValue), (uint8_t*) &psValue);
                    bam_aux_append(aln, "PQ", 'i', sizeof(pqValue), (uint8_t*) &pqValue);
                    totalTagCuonnt++;
                }
                else{
                    totalUnTagCuonnt++;
                }
            }
            result = sam_write1(out, bamHdr, aln);
            free(tmp.quality);
            free(tmp.qseq);
            free(tmp.op);
            free(tmp.ol);
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

int HaplotagProcess::judgeHaplotype(Alignment align, std::string chrName, double percentageThreshold, std::ofstream *tagResult, int &pqValue){

    int hp1Count = 0;
    int hp2Count = 0;
    //record variants on this read
    std::map<int,int> variants;

    // Skip variants that are to the left of this read
    while( firstVariantIter != currentVariants.end() && (*firstVariantIter).first < align.refStart )
        firstVariantIter++;
    
    if( firstVariantIter == currentVariants.end() )
        return 0;
     
    // position relative to reference
    int ref_pos = align.refStart;
    // position relative to read
    int query_pos = 0;
    // translation char* to string;
    std::string qseq = align.qseq;
    // translation char* to string;
        
    // set variant start for current alignment
    std::map<int, RefAlt>::iterator currentVariantIter = firstVariantIter;

    //std::map<int, RefAlt>::reverse_iterator last = currentVariants.rbegin();
        
    // reading cigar to detect snp on this read
    for(int i = 0; i < align.cigar_len ; i++ ){
        int cigar_op = align.op[i];
        int length   = align.ol[i];
 
        // iterator next variant
        while( currentVariantIter != currentVariants.end() && (*currentVariantIter).first < ref_pos ){
            currentVariantIter++;
        }
        
        // CIGAR operators: MIDNSHP=X correspond 012345678
        // 0: alignment match (can be a sequence match or mismatch)
        // 7: sequence match
        // 8: sequence mismatch
        if( cigar_op == 0 || cigar_op == 7 || cigar_op == 8 ){
                
            while( currentVariantIter != currentVariants.end() && (*currentVariantIter).first < ref_pos + length){
                
                int offset = (*currentVariantIter).first - ref_pos;
                
                if( offset < 0){
                }
                else{
                    std::string base = qseq.substr(query_pos + offset, 1);
                    
                    // base match SNP allele
                    if( base == (*currentVariantIter).second.Ref || base == (*currentVariantIter).second.Alt ){
                        
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
                                variants[(*currentVariantIter).first]=0;
                            }
                            if( base == chrVariantHP2[chrName][(*currentVariantIter).first]){
                                hp2Count++;
                                variants[(*currentVariantIter).first]=1;
                            }
                        }
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
            std::cerr<< "alignment find unsupported CIGAR operation from read: " << align.qname << "\n";
            exit(1);
        }
    }
    
    double min,max;
    
    auto readIter = readSVHapCount.find(align.qname);
    if( readIter != readSVHapCount.end() ){
        hp1Count += readSVHapCount[align.qname][0];
        hp2Count += readSVHapCount[align.qname][1];
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
        (*tagResult)<< align.qname       << "\t" 
                 << align.chr         << "\t" 
                 << align.refStart    << "\t" 
                 << max/(max+min)     << "\t" 
                 << hpResultStr       << "\t"
                 << hp1Count+hp2Count << "\t" 
                 << hp1Count          << "\t" 
                 << hp2Count          << "\t"
                 << pqValue           << "\t";
        // loop each variant         
        for(auto v : variants ){
            (*tagResult)<< " " << v.first << "," << v.second ;
        }
        (*tagResult)<< "\n";  
    }
    
    return hpResult;
}


HaplotagProcess::HaplotagProcess(HaplotagParameters params):
totalAlignment(0),totalSupplementary(0),totalSecondary(0),totalUnmapped(0),totalTagCuonnt(0),totalUnTagCuonnt(0),processBegin(time(NULL)),integerPS(false)
{
    std::cerr<< "phased SNP file:   " << params.snpFile             << "\n";
    std::cerr<< "phased SV file:    " << params.svFile              << "\n";
    std::cerr<< "input bam file:    " << params.bamFile             << "\n";
    std::cerr<< "output bam file:   " << params.resultPrefix+".bam" << "\n";
    std::cerr<< "number of threads: " << params.numThreads          << "\n";
    std::cerr<< "write log file:    " << (params.writeReadLog ? "true" : "false") << "\n";
    std::cerr<< "log file:          " << (params.writeReadLog ? (params.resultPrefix+".out") : "") << "\n";
    std::cerr<< "-------------------------------------------\n";
    std::cerr<< "filter mapping quality below:  " << params.qualityThreshold    << "\n";
    std::cerr<< "percentage threshold:          " << params.percentageThreshold << "\n";
    std::cerr<< "tag supplementary:             " << (params.tagSupplementary ? "true" : "false") << "\n";
    std::cerr<< "-------------------------------------------\n";

    // load SNP vcf file
    std::time_t begin = time(NULL);
    std::cerr<< "parsing SNP VCF ... ";
    parseSVFile = false;
    variantParser(params.snpFile);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // load SV vcf file
    if(params.svFile!=""){
        begin = time(NULL);
        std::cerr<< "parsing SV VCF ...\n";
        parseSVFile = true;
        variantParser(params.svFile);
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
    std::cerr<< "total tag alignment: " << totalTagCuonnt     << "\n";
    std::cerr<< "total untagged:      " << totalUnTagCuonnt   << "\n";
};
