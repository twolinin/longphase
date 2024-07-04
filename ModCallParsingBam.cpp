#include "ModCallParsingBam.h"
#include "PhasingGraph.h"
#include "Util.h"

#include <cmath>
#include <iostream>
#include <string.h>
#include <sstream>
#include <fstream>
#include <typeinfo>
#include "htslib/thread_pool.h"
#include "htslib/sam.h"

MethFastaParser::MethFastaParser(std::string fastaFile, std::vector<ReferenceChromosome> &chrInfo){
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
        chrInfo.emplace_back(seqname, faidx_fetch_seq(fai , seqname , 0 ,seqlen+1 , &ref_len), seqlen);
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

void MethBamParser::detectMeth(std::string chrName, int chr_len, int numThreads, std::vector<ReadVariant> &readVariantVec){
    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};

    for( auto bamFile: params->bamFileVec ){

        // open cram file
        samFile *fp_in = hts_open(bamFile.c_str(),"r");
        // set reference
        hts_set_fai_filename(fp_in, params->fastaFile.c_str());        
        // read header
        bam_hdr_t *bamHdr = sam_hdr_read(fp_in); 
        // initialize an alignment
        bam1_t *aln = bam_init1(); 
        hts_idx_t *idx = NULL;
        
        if ((idx = sam_index_load(fp_in, bamFile.c_str())) == 0) {
            std::cout<<"ERROR: Cannot open index for bam file\n";
            exit(1);
        }
        
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
                continue;
            }
            parse_CIGAR(*bamHdr,*aln, readVariantVec);

        }
        hts_idx_destroy(idx);
        bam_hdr_destroy(bamHdr);
        bam_destroy1(aln);
        sam_close(fp_in);
        hts_tpool_destroy(threadPool.pool);
    }
}

void MethBamParser::parse_CIGAR(const bam_hdr_t &bamHdr,const bam1_t &aln, std::vector<ReadVariant> &readVariantVec){
    // see detail https://github.com/samtools/htslib/blob/develop/sam_mods.c
    
    // hts_base_mod_state can parse MM, ML and MN tags
    hts_base_mod_state *mod_state = hts_base_mod_state_alloc();
    if( bam_parse_basemod( &aln, mod_state) < 0 ){
        std::cout<<bam_get_qname(&aln)<<"\tFail pase MM\\ML\n";
    }
    
    /* 
    // see detail https://github.com/samtools/htslib/issues/1550
    
    // Assuming there are a total of 10 types of modifications.
    n_mods = 10;
    hts_base_mod allmod[n_mods];
    while ((n = bam_next_basemod(&aln, mod_state, allmod, n_mods, &pos)) > 0) {
        // Report 'n'th mod at sequence position 'pos'
        for (j = 0; j < n && j < n_mods; j++) {
            //5mC  code:m ascii:109
            //5hmC code:h ascii:104
            //see sam tags pdf
            if (allmod[j].modified_base == 109) {
                queryMethVec.emplace_back(pos, allmod[j].qual);
            }
        }
    }*/
 
    ReadVariant *tmpReadResult = new ReadVariant();
    (*tmpReadResult).read_name = bam_get_qname(&aln);
    (*tmpReadResult).source_id = bamHdr.target_name[aln.core.tid];
    (*tmpReadResult).reference_start = aln.core.pos;
    (*tmpReadResult).is_reverse = bam_is_rev(&aln);
    
    int refstart = aln.core.pos;
    // forward strand is checked from head to tail
    // reverse strand is checked from tail to head
    int refpos = (bam_is_rev(&aln) ? refstart + 1 : refstart);
    //auto qmethiter = (align.is_reverse ? align.queryMethVec.end()-1 : align.queryMethVec.begin() );
    //std::cout<<(*tmpReadResult).read_name<<"\t"<<align.queryMethVec.size()<<"\n";
    int querypos = 0;
    
    hts_base_mod mods[5];
    int pos;
    // bam_next_basemod can iterate over this cached state.
    int n = bam_next_basemod(&aln, mod_state, mods, 5, &pos);
    if( n <= 0 ){ 
        hts_base_mod_state_free(mod_state);
        return;
    }
    
    //Parse CIGAR 
    for(int cigaridx = 0; cigaridx < int(aln.core.n_cigar) ; cigaridx++){
      uint32_t *cigar = bam_get_cigar(&aln);
      int cigar_op = bam_cigar_op(cigar[cigaridx]);
      int length = bam_cigar_oplen(cigar[cigaridx]);
      
      // CIGAR operators: MIDNSHP=X correspond 012345678
      // 0: alignment match (can be a sequence match or mismatch)
      // 7: sequence match
      // 8: sequence mismatch
      
      if(cigar_op == 0 || cigar_op == 7 || cigar_op == 8){
          while(true){
               if( pos > (querypos+length) ){
                    break;
                }
                if( n <= 0 ){
                    break;
                }
                int methrpos;
                if(bam_is_rev(&aln)){
                    methrpos = pos - querypos + refpos - 1;
                    
                }else{
                    methrpos = pos - querypos + refpos;
                }

                if( int((*refString).length()) < methrpos ){
                    break;
                }
                for (int j = 0; j < n && j < 5; j++) {
                    if (mods[j].modified_base == 109 && pos <= (querypos+length)) {
                        //modification
                        if( mods[j].qual >= params->modThreshold*255 ){
                            (*chrMethMap)[methrpos].methreadcnt++;
                            //strand - is 1, strand + is 0
                            (*chrMethMap)[methrpos].strand = (bam_is_rev(&aln) ? 1 : 0);
                            (*chrMethMap)[methrpos].modReadVec.push_back(bam_get_qname(&aln));
                            (*tmpReadResult).variantVec.emplace_back(methrpos, 0, 60);

                        }
                        //non-modification (not include no detected)
                        else if( mods[j].qual <= params->unModThreshold*255 ){ 
                            (*chrMethMap)[methrpos].canonreadcnt++;
                            (*chrMethMap)[methrpos].nonModReadVec.push_back(bam_get_qname(&aln));
                            (*tmpReadResult).variantVec.emplace_back(methrpos, 1, 60  );
                        }
                        else{
                            (*chrMethMap)[methrpos].noisereadcnt++;
                        }
                      }
                  }
                  n = bam_next_basemod(&aln, mod_state, mods, 5, &pos);
             }
          querypos += length;
          refpos += length;
      }
      // 1: insertion to the reference
      else if(cigar_op == 1){ 
          while(n>0 && pos <= (querypos+length)){
            n = bam_next_basemod(&aln, mod_state, mods, 5, &pos);
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
          while(n>0 && pos <= (querypos+length)){
            n = bam_next_basemod(&aln, mod_state, mods, 5, &pos);
          }
          querypos += length;
      }
      // 5: hard clipping (clipped sequences NOT present in SEQ)
      // 6: padding (silent deletion from padded reference)
      else if(cigar_op == 5 || cigar_op == 6){
          //do nothing
          refpos += length;
      }
    }
    
    hts_base_mod_state_free(mod_state);
    
    int refend = (bam_is_rev(&aln) ? refpos : refpos + 1);

    if(bam_is_rev(&aln)){
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

void MethBamParser::exportResult(std::string chrName, std::string chrSquence, int chrLen , std::map<int,std::vector<int>> &passPosition, std::ostringstream &modCallResult){
    
    for(std::map<int , MethPosInfo>::iterator posinfoIter = chrMethMap->begin(); posinfoIter != chrMethMap->end(); posinfoIter++){
        std::string infostr= "";
        std::string eachpos;
        std::string samplestr;
        std::string strandinfo;
        std::string ref;
        bool print = false;
        
        auto passPosIter =  passPosition.find((*posinfoIter).first);
        auto prepassPosIter =  passPosition.find((*posinfoIter).first - 1);
        auto nextpassPosIter =  passPosition.find((*posinfoIter).first + 1);
        
        if(!passPosition[passPosIter->first].empty() && !passPosition[prepassPosIter->first].empty()){
            for(std::vector<int>::iterator posIter = passPosition[passPosIter->first].begin(); posIter != passPosition[passPosIter->first].end(); posIter++){
                for(std::vector<int>::iterator posIter2 = passPosition[prepassPosIter->first].begin(); posIter2 != passPosition[prepassPosIter->first].end(); posIter2++){
                    if((*posIter)-1 == (*posIter2)){
                        print = true;
                    }
                }
            }
        }
        if(!passPosition[passPosIter->first].empty() && !passPosition[nextpassPosIter->first].empty()){
            for(std::vector<int>::iterator posIter = passPosition[passPosIter->first].begin(); posIter != passPosition[passPosIter->first].end(); posIter++){
                for(std::vector<int>::iterator posIter2 = passPosition[nextpassPosIter->first].begin(); posIter2 != passPosition[nextpassPosIter->first].end(); posIter2++){
                    if((*posIter)+1 == (*posIter2)){
                        print = true;
                    }
                }
            }
        }

        // prevent variant coordinates from exceeding the reference.
        if( chrLen < (*posinfoIter).first ){
            continue;
        }

        // prevent abnormal REF allele.
        ref = chrSquence.substr((*posinfoIter).first,1);
        if( ref != "A" && ref != "T" && ref != "C" && ref != "G" && 
            ref != "a" && ref != "t" && ref != "c" && ref != "g"  ){
            continue;
        }
        
        // set variant strand
        if( (*posinfoIter).second.strand == 1 ){
            strandinfo = "RS=N;";
        }
        else if( (*posinfoIter).second.strand == 0 ){
            strandinfo = "RS=P;";
        }
        else{
            continue;
        }
        
        //Output contains only consecutive methylation position (CpG)
        if( (prepassPosIter == passPosition.end() && nextpassPosIter == passPosition.end()) || passPosIter ==  passPosition.end() || print == false){
            if( params->outputAllMod ){
                int nonmethcnt = (*posinfoIter).second.canonreadcnt;
                samplestr = (*posinfoIter).second.heterstatus + ":" + std::to_string((*posinfoIter).second.methreadcnt) + ":" + std::to_string(nonmethcnt) + ":" + std::to_string((*posinfoIter).second.depth);

                eachpos = chrName + "\t" + std::to_string((*posinfoIter).first + 1) + "\t" + "." + "\t" + ref + "\t" + "N" + "\t" + "." + "\t" + "." + "\t" + strandinfo + infostr + "\t" + "GT:MD:UD:DP" + "\t" + samplestr + "\n";
                modCallResult<<eachpos;
            }         
            continue;
        }
        
        // append modification reads
        if((*posinfoIter).second.modReadVec.size() > 0 ){
            infostr += "MR=";
            for(auto readName : (*posinfoIter).second.modReadVec ){
                infostr += readName + ",";
            }
            infostr.back() = ';';
        }
            
        // append non modification reads
        if((*posinfoIter).second.nonModReadVec.size() > 0 ){
            infostr += "NR=";
            for(auto readName : (*posinfoIter).second.nonModReadVec ){
                infostr += readName + ",";
            }
            infostr.back() = ';';
        }

        if((*posinfoIter).second.heterstatus == "0/1"){
            int nonmethcnt = (*posinfoIter).second.canonreadcnt;
            samplestr = (*posinfoIter).second.heterstatus + ":" + std::to_string((*posinfoIter).second.methreadcnt) + ":" + std::to_string(nonmethcnt) + ":" + std::to_string((*posinfoIter).second.depth);

            eachpos = chrName + "\t" + std::to_string((*posinfoIter).first + 1) + "\t" + "." + "\t" + ref + "\t" + "N" + "\t" + "." + "\t" + "PASS" + "\t" + strandinfo + infostr + "\t" + "GT:MD:UD:DP" + "\t" + samplestr + "\n";

            modCallResult<<eachpos;
        }
    }
}

void writeResultVCF( ModCallParameters &params, std::vector<ReferenceChromosome> &chrInfo, std::map<std::string,std::ostringstream> &chrModCallResult){

    std::ofstream modCallResultVcf(params.resultPrefix+".vcf", std::ios_base::app);
    if(!modCallResultVcf.is_open()){
        std::cerr<<"Fail to open output file :\n";
    }
    else{
        // set vcf header
        modCallResultVcf<<"##fileformat=VCFv4.2\n";
        modCallResultVcf<<"##INFO=<ID=RS,Number=.,Type=String,Description=\"Read Strand\">\n";
        modCallResultVcf<<"##INFO=<ID=MR,Number=.,Type=String,Description=\"Read Name of Modified position\">\n";
        modCallResultVcf<<"##INFO=<ID=NR,Number=.,Type=String,Description=\"Read Name of nonModified position\">\n";
        modCallResultVcf<<"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        modCallResultVcf<<"##FORMAT=<ID=MD,Number=1,Type=Integer,Description=\"Modified Depth\">\n";
        modCallResultVcf<<"##FORMAT=<ID=UD,Number=1,Type=Integer,Description=\"Unmodified Depth\">\n";
        modCallResultVcf<<"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
        for(const auto& chrIter : chrInfo){
            modCallResultVcf<<"##contig=<ID="<<chrIter.name<<",length="<<chrIter.length<<">\n";
        }
        modCallResultVcf << "##longphaseVersion=" << params.version << "\n";
        modCallResultVcf << "##commandline=\""    << params.command << "\"\n";
        modCallResultVcf<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
        
        // set vcf body
        for(const auto& chrIter : chrInfo){
            modCallResultVcf << chrModCallResult[chrIter.name].str();
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
        }
        else if(methcnt >= nonmethcnt){
            (*chrmethmapIter).second.heterstatus = "1/1";
        }
        else{
            (*chrmethmapIter).second.heterstatus = "0/0";
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

MethylationGraph::MethylationGraph(ModCallParameters &in_params){
    params=&in_params;
    
    nodeInfo = new std::map<int,ReadBaseMap*>;
    edgeList = new std::map<int,VariantEdge*>;
    forwardModNode = new std::map<int,ReadBaseMap*>;
    reverseModNode = new std::map<int,ReadBaseMap*>;
}

MethylationGraph::~MethylationGraph(){
}

void MethylationGraph::addEdge(std::vector<ReadVariant> &in_readVariant){
    readVariant = &in_readVariant;
    // iter all read
    for(std::vector<ReadVariant>::iterator readIter = in_readVariant.begin() ; readIter != in_readVariant.end() ; readIter++ ){
        ReadVariant tmpRead;
        
        for( auto variant : (*readIter).variantVec ){
            // modification in the forward strand
            if( variant.quality == -2 ){
                auto nodeIter = forwardModNode->find(variant.position);
                if( nodeIter == forwardModNode->end() ){
                    (*forwardModNode)[variant.position] = new ReadBaseMap();
                }
                
                (*(*forwardModNode)[variant.position])[(*readIter).read_name] = 60;

                continue;
            }
            // modification in the reverse strand
            if( variant.quality == -3 ){
                auto nodeIter = reverseModNode->find(variant.position);
                if( nodeIter == reverseModNode->end() ){
                    (*reverseModNode)[variant.position] = new ReadBaseMap();
                }
                
                (*(*reverseModNode)[variant.position])[(*readIter).read_name] = 60;
                
                continue;
            }
            
            
            tmpRead.variantVec.push_back(variant);
            
            auto nodeIter = nodeInfo->find(variant.position);
            
            if( nodeIter == nodeInfo->end() ){
                (*nodeInfo)[variant.position] = new ReadBaseMap();
            }
            
            (*(*nodeInfo)[variant.position])[(*readIter).read_name] = variant.quality;
        }
        // iter all pair of snp and construct initial graph
        std::vector<Variant>::iterator variant1Iter = tmpRead.variantVec.begin();
        std::vector<Variant>::iterator variant2Iter = std::next(variant1Iter,1);
        while(variant1Iter != tmpRead.variantVec.end() && variant2Iter != tmpRead.variantVec.end() ){

            // create new edge if not exist
            std::map<int,VariantEdge*>::iterator posIter = edgeList->find((*variant1Iter).position);
            if( posIter == edgeList->end() )
                (*edgeList)[(*variant1Iter).position] = new VariantEdge((*variant1Iter).position);

            // add edge process
            for(int nextNode = 0 ; nextNode < params->connectAdjacent; nextNode++){

                // this allele support ref
                if( (*variant1Iter).allele == 0 )
                    (*edgeList)[(*variant1Iter).position]->ref->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).read_name,0,1);
                // this allele support alt
                if( (*variant1Iter).allele == 1 )
                    (*edgeList)[(*variant1Iter).position]->alt->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).read_name,0,1);
                
                // next snp
                variant2Iter++;
                if( variant2Iter == tmpRead.variantVec.end() ){
                    break;
                }
            }

            variant1Iter++;
            variant2Iter = std::next(variant1Iter,1);
        }
    }
}

void MethylationGraph::connectResults(std::string chrName, std::map<int,std::vector<int>> &passPosition){

    // check clear connect variant
    for(std::map<int,ReadBaseMap*>::iterator nodeIter = nodeInfo->begin() ; nodeIter != nodeInfo->end() ; nodeIter++ ){

        // check next position
        std::map<int,ReadBaseMap*>::iterator nextNodeIter = std::next(nodeIter, 1);
        if( nextNodeIter == nodeInfo->end() ){
             break;
        }

        int currPos = nodeIter->first;

        // Check if there is no edge from current node
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList->find( currPos );
        if( edgeIter==edgeList->end() )
            continue;

        // check connect between surrent SNP and next n SNPs
        for(int i = 0 ; i < params->connectAdjacent ; i++ ){
            int nextPos = nextNodeIter->first;
            // get number of RR read and RA read
            std::pair<int,int> tmp = edgeIter->second->findNumberOfRead(nextPos);
            int totalConnectReads = tmp.first + tmp.second;
            int minimumConnection = (((*(*nodeIter).second).size() + (*(*nextNodeIter).second).size())/4);

            double majorRatio = (double)std::max(tmp.first,tmp.second)/(double)(tmp.first+tmp.second);
            

            if( majorRatio >= params->connectConfidence && totalConnectReads > minimumConnection && tmp.first + tmp.second > 6 ){
                passPosition[currPos].push_back(nextPos);
                passPosition[nextPos].push_back(currPos);
            }
            nextNodeIter++;
            if( nextNodeIter == nodeInfo->end() )
                break;
        }
    }
}

void MethylationGraph::destroy(){

    for( auto edgeIter = edgeList->begin() ; edgeIter != edgeList->end() ; edgeIter++ ){
        edgeIter->second->ref->destroy();
        edgeIter->second->alt->destroy();
        delete edgeIter->second->ref;
        delete edgeIter->second->alt;
    }
    
    for( auto nodeIter = nodeInfo->begin() ; nodeIter != nodeInfo->end() ; nodeIter++ ){
        delete nodeIter->second;
    }
    
    for( auto nodeIter = forwardModNode->begin() ; nodeIter != forwardModNode->end() ; nodeIter++ ){
        delete nodeIter->second;
    }
    
    for( auto nodeIter = reverseModNode->begin() ; nodeIter != reverseModNode->end() ; nodeIter++ ){
        delete nodeIter->second;
    }
    
    delete nodeInfo;
    delete edgeList;
    delete forwardModNode;
    delete reverseModNode;
}
