#include "ModCallParsingBam.h"
#include "Util.h"

#include <cmath>
#include <iostream>
#include <string.h>
#include <sstream>
#include <fstream>
#include <typeinfo>
#include "htslib/thread_pool.h"
#include "htslib/sam.h"

MethFastaParser::MethFastaParser(std::string fastaFile, std::map<std::string, std::string> &chrString){
    if(fastaFile==""){
        return;
    }
    
    faidx_t *fai = NULL;
    fai = fai_load(fastaFile.c_str());
    int fai_nseq = faidx_nseq(fai);
    const char *seqname;
    int seqlen = 0;
    for(int i=0;i<fai_nseq;i++){
        int ref_len = 0;
        seqname = faidx_iseq(fai, i);
        seqlen = faidx_seq_len(fai, seqname);
        chrString[std::string(seqname)] = faidx_fetch_seq(fai , seqname , 0 ,seqlen+1 , &ref_len);
    }
}

MethFastaParser::~MethFastaParser(){
}

MethBamParser::MethBamParser(ModCallParameters &in_params, std::string &in_refString){
	params=&in_params;
	refString=&in_refString;
	refstartpos = 0;
    
    chrMethMap = new std::map<int , MethPosInfo>;
    readStartEndMap = new std::map<int, std::pair<int,int>>; 
}

MethBamParser::~MethBamParser(){
    delete chrMethMap;
    delete readStartEndMap;
}

static char *code(int id) {
    static char code[20];
    if (id > 0) {
        code[0] = id;
        code[1] = 0;
    }

    return code;
}
void MethBamParser::detectMeth(std::string chrName, int chr_len, std::vector<ReadVariant> &readVariantVec){
	int numThreads = params->numThreads;
    //char *bamList[numThreads];
    //int chunkSize = chr_len/numThreads + numThreads;
	
    // create list. i.e. chr1:100-200
    //for(int i=0;i<numThreads;i++){
    //    std::string tmp = chrName + ":" + std::to_string(i*chunkSize+1) + "-" + std::to_string((i+1)*chunkSize);
    //    bamList[i] = new char[tmp.length()+1];
    //    strcpy(bamList[i], tmp.c_str());
    //}
    
    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};

	// open cram file
	samFile *fp_in = hts_open(params->methylBamFile.c_str(),"r");
	// set reference
	hts_set_fai_filename(fp_in, params->fastaFile.c_str());		
    // read header
	bam_hdr_t *bamHdr = sam_hdr_read(fp_in); 
    // initialize an alignment
	bam1_t *aln = bam_init1(); 
    hts_idx_t *idx = NULL;
    
    if ((idx = sam_index_load(fp_in, params->methylBamFile.c_str())) == 0) {
        std::cout<<"ERROR: Cannot open index for bam file\n";
        exit(1);
    }
    
    /*
    hts_itr_multi_t *iter;
    if( (iter = sam_itr_regarray(idx, bamHdr, bamList, numThreads)) == 0){
        std::cout<<"ERROR: Cannot open iterator for region %s for bam file\n";
        exit(1);
    }
    */
    std::string range = chrName + ":1-" + std::to_string(chr_len);
    hts_itr_t* iter = sam_itr_querys(idx, bamHdr, range.c_str());

    
    int result;
    
    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
    }
    hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &threadPool);

    while ((result = sam_itr_multi_next(fp_in, iter, aln)) >= 0) { 
        int flag = aln->core.flag;

        if (  aln->core.qual < 1  // mapping quality
             || (flag & 0x4)   != 0  // read unmapped
             || (flag & 0x100) != 0  // secondary alignment. repeat. 
                                     // A secondary alignment occurs when a given read could align reasonably well to more than one place.
             || (flag & 0x400) != 0  // duplicate 
             || (flag & 0x800) != 0  // supplementary alignment
                                     // A chimeric alignment is represented as a set of linear alignments that do not have large overlaps.
             ){ 
            std::string str = bamHdr->target_name[aln->core.tid];
            continue;
        }

        AlignmentMethExtend *tmp = new AlignmentMethExtend;

        tmp->chr = bamHdr->target_name[aln->core.tid];
        tmp->refStart = aln->core.pos;
        tmp->qname = bam_get_qname(aln);
        tmp->qlen = aln->core.l_qseq;
        tmp->cigar_len = aln->core.n_cigar;
        tmp->op = (int*)malloc(tmp->cigar_len*sizeof(int));
        tmp->ol = (int*)malloc(tmp->cigar_len*sizeof(int));
        tmp->is_reverse = bam_is_rev(aln);
        
        memset(tmp->op, 0, tmp->cigar_len);
        memset(tmp->ol, 0, tmp->cigar_len);
        
        // set string size
        tmp->qseq = (char *)malloc( aln->core.l_qseq + 1 );
        memset(tmp->qseq, 0, tmp->qlen+1);
        // set string size
        tmp->quality = (char *)malloc( aln->core.l_qseq + 1 );
        memset(tmp->quality, 0, tmp->qlen+1);
        
        //quality string
        uint8_t *q = bam_get_seq(aln); 
        uint8_t *quality = bam_get_qual(aln);
        
        for(int i=0; i< aln->core.l_qseq ; i++){
			// gets nucleotide id and converts them into IUPAC id.
            tmp->qseq[i] = seq_nt16_str[bam_seqi(q,i)]; 
            // get base quality
            tmp->quality[i] = quality[i];
		}
        
        uint32_t *cigar = bam_get_cigar(aln);
        // store cigar
        for(unsigned int k =0 ; k < aln->core.n_cigar ; k++){
            tmp->op[k] = bam_cigar_op(cigar[k]);
            tmp->ol[k] = bam_cigar_oplen(cigar[k]);
        }
		
        //get "ML" and "MM" tag
        tmp->mmStr = get_aux_tag(aln, "MM");
        tmp->mlStr = get_aux_tag(aln, "ML");
		//std::cout<<tmp->qname<<"\t"<<tmp->mmStr<<"\t"<<tmp->mlStr<<"\n";
		
		
		hts_base_mod_state *m = hts_base_mod_state_alloc();
		if(bam_parse_basemod(aln, m)<0){
			std::cout<<tmp->qname<<"\tFail pase MM\ML\n";
		}
        
		int k,j,n,pos;
		hts_base_mod mods[5];
		
		while ((n=bam_next_basemod(aln, m, mods, 5, &pos)) > 0) {
            
            for (j = 0; j < n && j < 5; j++) {
                int m_strand, m_implicit;
				char m_canonical;
				
				int ret = bam_mods_query_type(m, mods[j].modified_base,
                                                  &m_strand, &m_implicit,
                                                  &m_canonical);
				/*std::cout<<tmp->qname<<"\t"<<pos<<"\t"<<
											seq_nt16_str[bam_seqi(bam_get_seq(aln), pos)]<<"\t"<<
											mods[j].canonical_base<<"\t"<<
											mods[j].strand<<"\t"<<
											code(mods[j].modified_base)<<"\t"<<
											"?."[m_implicit]<<"\t"<<
											mods[j].qual<<"\n";*/
											
				if(mods[j].modified_base==109 ){
					MethPosProb *mpp = new MethPosProb(pos, mods[j].qual);
					tmp->queryMethVec.push_back((*mpp)); 
					delete mpp;
				}
            }
			
		}
		
		hts_base_mod_state_free(m);
        /*if( tmp->mmStr.length() != 0 && tmp->mlStr.length() != 0 ){
            getmeth((*tmp));
			parse_CIGAR((*tmp), readVariantVec);
        }*/
		
		parse_CIGAR((*tmp), readVariantVec);
        //release memory
        tmp->queryMethVec.clear();
        tmp->queryMethVec.shrink_to_fit();
        free(tmp->quality);
        free(tmp->qseq);
        free(tmp->op);
        free(tmp->ol);
        delete tmp;
	}
	hts_idx_destroy(idx);
    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
	sam_close(fp_in);
    hts_tpool_destroy(threadPool.pool);
}

std::string MethBamParser::get_aux_tag(const bam1_t *aln, const char tag[2]){
	kstring_t res = KS_INITIALIZE;
    //if there have "MM" tag, then convert to string type
    if(bam_aux_get_str(aln, tag, &res) == 1){ 
        int len = ks_len(&res);
        char *ks_s = ks_str(&res);
        //because kstring don't have '\0' stop charcter
        std::string s(len, '\0'); 
        for(int i=0;i<len;i++){
            s[i] = ks_s[i];
        }
        ks_free(&res);
        return s;
    }
    else{
        ks_free(&res);
        return "";
    }
}

void MethBamParser::getmeth(AlignmentMethExtend &align){
	
	//transalte MM and ML TAG to the reference genome
	//eg: MM: 0, 1, 3  ; ML: 23, 0, 250
	//eg: SEQ: AGCTTTCAGCCTAGCCCTTTCA ; strand = +
	//it means AG'C'TTTCAG'C'CTAGCC'C'TTTCA ==> left to right, each 'C' probability is : 23/255, 0/255, 250/255
	//eg: MM: 0, 1, 1  ; ML: 67, 100, 200
	//eg: SEQ: CGATTTAGGACATAGAAAAGAT ; strand = -
	//it means C'G'ATTTAG'G'ACATAGAAAA'G'AT ==> right to left, each 'G' probability is : 67/255, 100/255, 200/255
	
	//parse mmtag (methylation position)
    //position
    std::vector<int> mmVec; 
    //create string stream from the string
    std::stringstream mm_stream(align.mmStr); 
    while(mm_stream.good()) {
        std::string substr;
        //get first string delimited by comma
        getline(mm_stream, substr, ','); 
		if(substr.find(";C+m")!=std::string::npos){
			//mmVec.push_back(std::stoi(substr));
			break;
		}
        if(substr.find("MM") == std::string::npos){
            mmVec.push_back(std::stoi(substr));
        }
		
    }
	
	//parse mltag (methylation probability)
    //prob
    std::vector<float> mlVec; 
    //create string stream from the string
    std::stringstream ml_stream(align.mlStr); 
    while(ml_stream.good()) {
        std::string substr;
        //get first string delimited by comma
        getline(ml_stream, substr, ','); 
        if(substr.find("ML") == std::string::npos){
            mlVec.push_back(std::stof(substr));
        }
    }
	
    if(int(mmVec.size()) == 0 || int(mlVec.size()) == 0){
		//std::cerr<<"There is no MM Tag or ML Tag\n";
        return;
    }

    std::vector<int>::iterator mmiter = mmVec.begin();
    std::vector<float>::iterator mliter = mlVec.begin();
	for(int i=0;i<mlVec.size()/2;i++){
		mliter++;
	}
	//std::cout<<align.qname<<"\t"<<mmVec.size()<<"\t"<<mlVec.size()<<"\n";
	if(!align.is_reverse){
		int seqi = 0;
		for( ; mmiter != mmVec.end(), mliter != mlVec.end(); mmiter++, mliter++){
			if(mmiter == mmVec.end()){
				break;
			}
			int intervalC = (*mmiter) + 1;
			while(seqi < align.qlen && intervalC > 0){
				if(align.qseq[seqi] == 'C'){
					intervalC--;
				}
				seqi++;
			}
			MethPosProb *mpp = new MethPosProb(seqi, *mliter);
			align.queryMethVec.push_back((*mpp));
            delete mpp;
		}
	}
	else{
		int seqi_rev = align.qlen - 1;
		for( ; mmiter != mmVec.end(), mliter != mlVec.end(); mmiter++, mliter++){
			if(mmiter == mmVec.end()){
				break;
			}
			int intervalC = (*mmiter) + 1;
			while(seqi_rev >= 0 && intervalC > 0){
				if(align.qseq[seqi_rev] == 'G'){
					intervalC--;
				}
				seqi_rev--;
			}
			MethPosProb *mpp = new MethPosProb(seqi_rev+1, *mliter);
			align.queryMethVec.push_back((*mpp)); 
            delete mpp;
		}
	}
    
    mmVec.clear();
    mmVec.shrink_to_fit();
	mlVec.clear();
    mlVec.shrink_to_fit();
}

void MethBamParser::parse_CIGAR(AlignmentMethExtend &align, std::vector<ReadVariant> &readVariantVec){

	if(align.queryMethVec.size() == 0){
		return;
	}
    //---------------------------------------------------------------------
    ReadVariant *tmpReadResult = new ReadVariant();
    (*tmpReadResult).read_name = align.qname;
    (*tmpReadResult).source_id = align.chr;
    (*tmpReadResult).reference_start = align.refStart;
    (*tmpReadResult).is_reverse = align.is_reverse;
    //---------------------------------------------------------------------
    
    
    int refstart = align.refStart;
	int refend = refstart + 1;
    // forward strand is checked from head to tail
    // reverse strand is checked from tail to head
    int refpos = (align.is_reverse ? refstart + 1 : refstart);
    //auto qmethiter = (align.is_reverse ? align.queryMethVec.end()-1 : align.queryMethVec.begin() );
	auto qmethiter = align.queryMethVec.begin() ;
    //std::cout<<(*tmpReadResult).read_name<<"\t"<<align.queryMethVec.size()<<"\n";
    int querypos = 0;
    //Parse CIGAR 
    for(int cigaridx = 0; cigaridx < align.cigar_len ; cigaridx++){
        
        int cigar_op = align.op[cigaridx];
        int length = align.ol[cigaridx];
        
        // CIGAR operators: MIDNSHP=X correspond 012345678
        // 0: alignment match (can be a sequence match or mismatch)
        // 7: sequence match
        // 8: sequence mismatch
        if(cigar_op == 0 || cigar_op == 7 || cigar_op == 8){
            
            while(true){
                
                if( (*qmethiter).position > (querypos+length) ){
                    break;
                }
                else if( qmethiter == (align.is_reverse ? align.queryMethVec.begin() : align.queryMethVec.end()) ){
                    break;
                }
                //translate query (read) position to reference position
                int methrpos = (*qmethiter).position - querypos + refpos;
				//std::cout<<(*refString).length()<<"\t"<<methrpos<<"\n";
				if((*refString).length()<methrpos){
					break;
				}
				std::string refChar = (*refString).substr(methrpos-1,1);
				//std::cout<<"parser_CIGA\t"<<align.qname<<"\t"<<methrpos<<"\t"<<(*qmethiter).prob<<"\t"<<(*qmethiter).position<<"\t"<<querypos<<"\t"<<refpos<<"\n";
                if( (refChar == "C") ? !align.is_reverse : align.is_reverse ){
					//methylation
					if( (*qmethiter).prob >= params->modThreshold*255 ){
						(*chrMethMap)[methrpos].methreadcnt++;
						//strand - is 1, strand + is 0
						(*chrMethMap)[methrpos].strand = (align.is_reverse ? 1 : 0); 
						(*chrMethMap)[methrpos].modReadVec.push_back(align.qname);
						//---------------------------------------------------------------------
						Variant *tmpVariant = new Variant(methrpos, 0, 60 );
						(*tmpReadResult).variantVec.push_back( (*tmpVariant) );
						delete tmpVariant;
						//---------------------------------------------------------------------
					}
					//non-methylation (not include no detected)
					else if( (*qmethiter).prob <= params->unModThreshold*255 ){ 
						(*chrMethMap)[methrpos].canonreadcnt++;
						(*chrMethMap)[methrpos].nonModReadVec.push_back(align.qname);
						//---------------------------------------------------------------------
						Variant *tmpVariant = new Variant(methrpos, 1, 60 );
						(*tmpReadResult).variantVec.push_back( (*tmpVariant) );
						delete tmpVariant;
						//---------------------------------------------------------------------
					}
					else{
						(*chrMethMap)[methrpos].noisereadcnt++;
					}	
				}
                if(align.is_reverse){
                    qmethiter--;
                }
                else{
                    qmethiter++;
                }
            }
            querypos += length;
            refpos += length;
        }
        // 1: insertion to the reference
        else if(cigar_op == 1){ 
            if(align.is_reverse){
                while(qmethiter != align.queryMethVec.begin() && (*qmethiter).position <= (querypos+length)){
                    qmethiter--;
                }
            }
            else{
                while(qmethiter != align.queryMethVec.end() && (*qmethiter).position <= (querypos+length)){
                    qmethiter++;
                }
            }
            querypos += length;
            
        }
        // 2: deletion from the reference
        else if(cigar_op == 2){
            refpos += length;
        }
        // 3: skipped region from the reference
        else if(cigar_op == 3){
            refpos += length;
        }
        // 4: soft clipping (clipped sequences present in SEQ)
        else if(cigar_op == 4){
            if(align.is_reverse){
                while(qmethiter != align.queryMethVec.begin() && (*qmethiter).position <= (querypos+length)){
                    qmethiter--;
                }
            }
            else{
                while(qmethiter != align.queryMethVec.end() && (*qmethiter).position <= (querypos+length)){
                    qmethiter++;
                }
            }
            querypos += length;
            
        }
        // 5: hard clipping (clipped sequences NOT present in SEQ)
        // 6: padding (silent deletion from padded reference)
        else if(cigar_op == 5 || cigar_op == 6){
            //do nothing
        }
    }
    
    refend = (align.is_reverse ? refpos : refpos + 1);

    if(align.is_reverse){
        (*readStartEndMap)[refstart+1].second += 1;
		(*readStartEndMap)[refend].second -= 1;
    }
    else{
        (*readStartEndMap)[refstart+1].first += 1;
        (*readStartEndMap)[refend].first -= 1;
    }
    
    if( (*tmpReadResult).variantVec.size() > 0 )
        readVariantVec.push_back((*tmpReadResult));
    
    delete tmpReadResult;
}

void MethBamParser::writeResultVCF( std::string chrName, std::map<std::string, std::string> &chrString, bool isFirstChr, std::map<int,int> &passPosition){

	std::ofstream ModResultVcf(params->resultPrefix+".vcf", std::ios_base::app);
	if(!ModResultVcf.is_open()){
		std::cerr<<"Fail to open output file :\n";
	}
    else{
		// set vcf header
        if(isFirstChr){ 
			ModResultVcf<<"##fileformat=VCFv4.2\n";
			ModResultVcf<<"##INFO=<ID=RS,Number=.,Type=String,Description=\"Read Strand\">\n";
			ModResultVcf<<"##INFO=<ID=MR,Number=.,Type=String,Description=\"Read Name of Modified position\">\n";
            ModResultVcf<<"##INFO=<ID=NR,Number=.,Type=String,Description=\"Read Name of nonModified position\">\n";
			ModResultVcf<<"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
			ModResultVcf<<"##FORMAT=<ID=MD,Number=1,Type=Integer,Description=\"Modified Depth\">\n";
			ModResultVcf<<"##FORMAT=<ID=UD,Number=1,Type=Integer,Description=\"Unmodified Depth\">\n";
			ModResultVcf<<"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\"\n";
			for(std::map<std::string, std::string>::iterator chrStringIter = chrString.begin(); chrStringIter != chrString.end(); chrStringIter++){
				ModResultVcf<<"##contig=<ID="<<(*chrStringIter).first<<",length="<<(*chrStringIter).second.length()<<">\n";
			}
			ModResultVcf<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
		}

		for(std::map<int , MethPosInfo>::iterator posinfoIter = chrMethMap->begin(); posinfoIter != chrMethMap->end(); posinfoIter++){
			std::string infostr= "";
            std::string eachpos;
            std::string samplestr;
            std::string strandinfo;
            
            int confidentModCount = 0;
            int confidentnonModCount = 0;
            
            auto passPosIter =  passPosition.find((*posinfoIter).first);
            
            if( passPosIter ==  passPosition.end() ){
                continue;
            }
            
			if((*posinfoIter).second.strand == -1){
				continue;
			}
            
            // append modification reads
            if((*posinfoIter).second.modReadVec.size() > 0 ){
                infostr += "MR=";
                for(auto readName : (*posinfoIter).second.modReadVec ){
                    confidentModCount++;
                    infostr += readName + ",";
                }
                infostr.back() = ';';
            }
                
            // append non modification reads
			if((*posinfoIter).second.nonModReadVec.size() > 0 ){
                infostr += "NR=";
                for(auto readName : (*posinfoIter).second.nonModReadVec ){
                    confidentnonModCount++;
                    infostr += readName + ",";
                }
                infostr.back() = ';';
            }
			
            if( confidentModCount == 0 || confidentnonModCount == 0 ){
                continue;
            }
            
			if((*posinfoIter).second.strand==1){
				strandinfo = "RS=N;";
			}
			else if((*posinfoIter).second.strand==0){
				strandinfo = "RS=P;";
			}
			
			if((*posinfoIter).second.heterstatus == "0/1"){
				int nonmethcnt = (*posinfoIter).second.canonreadcnt;
				samplestr = (*posinfoIter).second.heterstatus + ":" + std::to_string((*posinfoIter).second.methreadcnt) + ":" + std::to_string(nonmethcnt) + ":" + std::to_string((*posinfoIter).second.depth);
				if(chrString[chrName].length()<(*posinfoIter).first){
					break;
				}
				std::string ref = chrString[chrName].substr((*posinfoIter).first-1,1);
                eachpos = chrName + "\t" + std::to_string((*posinfoIter).first) + "\t" + "." + "\t" + ref + "\t" + "." + "\t" + "." + "\t" + "." + "\t" + strandinfo + infostr + "\t" + "GT:MD:UD:DP" + "\t" + samplestr + "\n";
				ModResultVcf<<eachpos;
			}
		}
	}
}

void MethBamParser::judgeMethGenotype(std::string chrName, std::vector<ReadVariant> &readVariantVec, std::vector<ReadVariant> &fReadVariantVec, std::vector<ReadVariant> &rReadVariantVec){
	
    for(std::map<int, MethPosInfo>::iterator chrmethmapIter = chrMethMap->begin(); chrmethmapIter != chrMethMap->end(); chrmethmapIter++){
		// methylation and noise are the feature to identify ASM candidate
		float methcnt      = (*chrmethmapIter).second.methreadcnt;
		float nonmethcnt   = (*chrmethmapIter).second.canonreadcnt;
        float depth        = (*chrmethmapIter).second.depth;
		float noisereadcnt = depth - methcnt - nonmethcnt;
        
        
		if(methcnt < 0 || nonmethcnt < 0){
			continue;
		}
		if(std::max(methcnt, nonmethcnt) == 0){
			continue;
		}

		float heterRatio = std::min(methcnt, nonmethcnt) / std::max(methcnt, nonmethcnt);
		float noiseRatio = noisereadcnt / depth;

		if(heterRatio >= params->heterRatio && noiseRatio <= params->noiseRatio ){
            (*chrmethmapIter).second.heterstatus = "0/1";
			 std::cout<< "heterstatus0/1\t" << chrName << "\t" <<(*chrmethmapIter).first << "\t" << heterRatio << "\t" << noiseRatio <<"\n";
        }
		else if(methcnt >= nonmethcnt){
			(*chrmethmapIter).second.heterstatus = "1/1";
			std::cout<< "heterstatus1/1\t" << chrName << "\t" <<(*chrmethmapIter).first << "\t" << heterRatio << "\t" << noiseRatio <<"\n";
		}
		else{
			(*chrmethmapIter).second.heterstatus = "0/0";
			std::cout<< "heterstatus0/0\t" << chrName << "\t" <<(*chrmethmapIter).first << "\t" << heterRatio << "\t" << noiseRatio <<"\n";
		}

	}
    
    for(std::vector<ReadVariant>::iterator rIter = readVariantVec.begin() ; rIter != readVariantVec.end() ; rIter++ ){

        ReadVariant fVec;
        ReadVariant rVec;
        
        for(std::vector<Variant>::iterator vIter = (*rIter).variantVec.begin(); vIter != (*rIter).variantVec.end() ; vIter++){
            if( (*chrMethMap)[(*vIter).position].heterstatus == "0/1" ){
                if((*rIter).is_reverse){
                    rVec.variantVec.push_back((*vIter));
                }
                else{
                    fVec.variantVec.push_back((*vIter));
                }
                //std::cout<< (*vIter).position << "\t" << (*rIter).is_reverse << "\n"; 
            }
        }
        
        if( fVec.variantVec.size() > 0 ){
            fVec.read_name = (*rIter).read_name;
            fVec.is_reverse = false;
            fReadVariantVec.push_back(fVec);
            
            fVec.variantVec.clear();
            fVec.variantVec.shrink_to_fit();
        }
        if( rVec.variantVec.size() > 0 ){
            rVec.read_name = (*rIter).read_name;
            rVec.is_reverse = true;
            rReadVariantVec.push_back(rVec);
            
            rVec.variantVec.clear();
            rVec.variantVec.shrink_to_fit();
        }
    }
    
}

void MethBamParser::calculateDepth(){
    std::map<int , MethPosInfo>::iterator methIter = chrMethMap->begin();
	std::pair<int,int> currdepth = std::make_pair(0,0);
	
	for(auto ReadIter = readStartEndMap->begin(); ReadIter != readStartEndMap->end(); ReadIter++){
		
		std::map<int, std::pair<int,int>>::iterator nextReadIter = std::next(ReadIter, 1);
		if(methIter == chrMethMap->end()){
			break;
		}
		if(nextReadIter == readStartEndMap->end()){
			break;
		}
		currdepth.first += (*ReadIter).second.first;
		currdepth.second += (*ReadIter).second.second;
		
		while((*methIter).first >= (*ReadIter).first && (*methIter).first < (*nextReadIter).first && methIter != chrMethMap->end()){
			//std::cout<<  (*methIter).first << "\t" << (*ReadIter).first << "\t" << (*nextReadIter).first << "\n";
            // strand +
            if((*methIter).second.strand == 0){ 
				(*methIter).second.depth = currdepth.first;
			}
            // strand -
			else if((*methIter).second.strand == 1){
				(*methIter).second.depth = currdepth.second;
			}
			
			methIter++;
		}
	}
	
	readStartEndMap->clear();
}
