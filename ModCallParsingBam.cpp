#include "ModCallParsingBam.h"
#include "PhasingGraph.h"
#include "Util.h"

#include <cmath>
#include <iostream>
#include <string.h>
#include <sstream>
#include <fstream>
#include <typeinfo>
#include <algorithm>
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

MethBamParser::MethBamParser(std::string inputChrName, ModCallParameters &in_params, MethSnpParser &snpMap, std::string &ref_string):chrName(inputChrName){
    params=&in_params;
    refString = &ref_string;
    refstartpos = 0;    
    chrMethMap = new std::map<int , MethPosInfo>;
    readStartEndMap = new std::map<int, std::pair<int,int>>; 

    currentVariants = new std::map<int, RefAlt>;

    if(snpMap.hasValidSnpData()){
        (*currentVariants) = snpMap.getVariants_markindel(chrName, ref_string);
    }

    firstVariantIter = currentVariants->begin();
}

MethBamParser::~MethBamParser(){
    delete chrMethMap;
    delete readStartEndMap;
    delete currentVariants;
}

void MethBamParser::detectMeth(std::string chrName, int lastSNPPos, htsThreadPool &threadPool, std::vector<ReadVariant> &readVariantVec){
    // record SNP start iter
    std::map<int, RefAlt>::iterator tmpFirstVariantIter = firstVariantIter;
    for( auto bamFile: params->bamFileVec ){
        firstVariantIter = tmpFirstVariantIter;
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
        
        std::string range = chrName + ":1-" + std::to_string(lastSNPPos);
        hts_itr_t* iter = sam_itr_querys(idx, bamHdr, range.c_str());

        
        int result;
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
    }
}

void MethBamParser::parse_CIGAR(const bam_hdr_t &bamHdr,const bam1_t &aln, std::vector<ReadVariant> &readVariantVec){
    // see detail https://github.com/samtools/htslib/blob/develop/sam_mods.c
    
    // hts_base_mod_state can parse MM, ML and MN tags
    hts_base_mod_state *mod_state = hts_base_mod_state_alloc();
    if( bam_parse_basemod( &aln, mod_state) < 0 ){
        hts_base_mod_state_free(mod_state);
        return;
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
    int ref_pos = aln.core.pos;
    //auto qmethiter = (align.is_reverse ? align.queryMethVec.end()-1 : align.queryMethVec.begin() );
    int querypos = 0;
    int aln_core_n_cigar = aln.core.n_cigar;
    uint32_t *cigar = bam_get_cigar(&aln);
    
    hts_base_mod mods[5];
    int pos;
    // bam_next_basemod can iterate over this cached state.
    int n = bam_next_basemod(&aln, mod_state, mods, 5, &pos);
    if( n <= 0 ){ 
        hts_base_mod_state_free(mod_state);
        delete tmpReadResult;
        return;
    }

    while( firstVariantIter != currentVariants->end() && (*firstVariantIter).first < ref_pos ){
        firstVariantIter++;
    }
    // set variant start for current alignment
    std::map<int, RefAlt>::iterator currentVariantIter = firstVariantIter;
    
    //Parse CIGAR 
    for(int cigaridx = 0; cigaridx < int(aln.core.n_cigar) ; cigaridx++){

        int cigar_op = bam_cigar_op(cigar[cigaridx]);
        int length = bam_cigar_oplen(cigar[cigaridx]);
        int variantPos = (*currentVariantIter).first;
        // get the first variant detected by the alignment.
        /*while( currentVariantIter != currentVariants->end() && variantPos < ref_pos ){
            currentVariantIter++;
            variantPos = (*currentVariantIter).first;
        }*/
        //while((currentVariantIter != currentVariants->end() && variantPos < ref_pos + length)){
            if( cigar_op == 0 || cigar_op == 7 || cigar_op == 8 ){
                while(currentVariantIter != currentVariants->end() && currentVariantIter->first < ref_pos + length){
                    int variantPos = currentVariantIter->first;
                    if(variantPos >= ref_pos){
                        int refAlleleLen = (*currentVariantIter).second.Ref.length();
                        int altAlleleLen = (*currentVariantIter).second.Alt.length();
                        int offset = variantPos - ref_pos;
                        int base_q = 0;
                        int allele = -1;

                        // The position of the variant exceeds the length of the read.
                        if( querypos + offset + 1 > int(aln.core.l_qseq) ){
                            return;
                        }

                        // SNP
                        if( refAlleleLen == 1 && altAlleleLen == 1){
                            char base = seq_nt16_str[bam_seqi(bam_get_seq(&aln), querypos + offset)];
                            if(base == (*currentVariantIter).second.Ref[0])
                                allele = 0;
                            else if(base == (*currentVariantIter).second.Alt[0])
                                allele = 1;

                            base_q = bam_get_qual(&aln)[querypos + offset];
                        } 
                
                        // insertion
                        //if( refAlleleLen == 1 && altAlleleLen != 1 && align.op[i+1] == 1 && i+1 < align.cigar_len){
                        if( refAlleleLen == 1 && altAlleleLen != 1 && cigaridx+1 < aln_core_n_cigar){
                    
                            // currently, qseq conversion is not performed. Below is the old method for obtaining insertion sequence.
                    
                            // uint8_t *qstring = bam_get_seq(aln); 
                            // qseq[i] = seq_nt16_str[bam_seqi(qstring,i)]; 
                            // std::string prevIns = ( align.op[i-1] == 1 ? qseq.substr(prev_query_pos, align.ol[i-1]) : "" );

                            if ( ref_pos + length - 1 == variantPos && bam_cigar_op(cigar[cigaridx+1]) == 1 ) {
                                allele = 1 ;
                            }
                            else {
                                allele = 0 ;
                            }
                            // using this quality to identify indel
                            base_q = -4;

                            // using this quality to identify danger indel
                            if ( (*currentVariantIter).second.is_danger ) {
                                base_q = -5 ;
                            }
                        } 
                
                        // deletion
                        //if( refAlleleLen != 1 && altAlleleLen == 1 && align.op[i+1] == 2 && i+1 < align.cigar_len){
                        if( refAlleleLen != 1 && altAlleleLen == 1 && cigaridx+1 < aln_core_n_cigar) {

                            if ( ref_pos + length - 1 == variantPos && bam_cigar_op(cigar[cigaridx+1]) == 2 ) {
                                allele = 1 ;
                            }
                            else {
                                allele = 0 ;
                            }
                            // using this quality to identify indel
                            base_q = -4;

                            // using this quality to identify danger indel
                            if ( (*currentVariantIter).second.is_danger ) {
                                base_q = -5 ;
                            }
                        } 
                
                        if( allele != -1){
                            // record snp result
                            (*tmpReadResult).variantVec.emplace_back(variantPos, allele, base_q, VariantType::SNP); 
                            (*chrMethMap)[variantPos].variantType = VariantType::SNP;                      
                        }
                    }
                    currentVariantIter++;
                }
            }
            //else break;
        //}

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
                    
                }
                else{
                    methrpos = pos - querypos + refpos;
                }

                if( int((*refString).length()) < methrpos ){
                    break;
                }

                for (int j = 0; j < n && j < 5; j++) {
                    auto it = chrMethMap->find(methrpos);
                    if (mods[j].modified_base == 109 && pos <= (querypos+length) && (it == chrMethMap->end() || it->second.variantType == VariantType::MOD)){
                        //modification
                        if( mods[j].qual >= params->modThreshold*255){
                            (*chrMethMap)[methrpos].methreadcnt++;
                            (*chrMethMap)[methrpos].variantType = VariantType::MOD;
                            //strand - is 1, strand + is 0
                            (*chrMethMap)[methrpos].strand = (bam_is_rev(&aln) ? 1 : 0);
                            (*chrMethMap)[methrpos].modReadVec.emplace_back(bam_get_qname(&aln));
                            (*tmpReadResult).variantVec.emplace_back(methrpos, 0, 60, VariantType::MOD);
                        }
                        //non-modification (not include no detected)
                        else if( mods[j].qual <= params->unModThreshold*255 ){ 
                            (*chrMethMap)[methrpos].canonreadcnt++;
                            (*chrMethMap)[methrpos].nonModReadVec.emplace_back(bam_get_qname(&aln));
                            (*tmpReadResult).variantVec.emplace_back(methrpos, 1, 60, VariantType::MOD);
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
            ref_pos += length;
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
            if(*refString != ""){
                int del_len = length;
                if ( ref_pos + del_len + 1 == (*currentVariantIter).first ){
                    //if( homopolymerLength((*currentVariantIter).first , ref_string) >=3 ){
                        // special case
                    //}
                }
                else if( (*currentVariantIter).first >= ref_pos  && (*currentVariantIter).first < ref_pos + del_len ){
                    // check snp in homopolymer
                    if( homopolymerLength((*currentVariantIter).first , *refString) >=3 ){
                        
                        int refAlleleLen = (*currentVariantIter).second.Ref.length();
                        int altAlleleLen = (*currentVariantIter).second.Alt.length();
                        int base_q = 0;  
                        
                        if( querypos + 1 > aln.core.l_qseq ){
                            hts_base_mod_state_free(mod_state);
                            delete tmpReadResult;
                            return;
                        }

                        int allele = -1;
                        // SNP
                        if( refAlleleLen == 1 && altAlleleLen == 1){
                            // get the next match
                            char base = seq_nt16_str[bam_seqi(bam_get_seq(&aln), querypos)];
                            if( base == (*currentVariantIter).second.Ref[0] ){
                                allele = 0;
                            }
                            else if( base == (*currentVariantIter).second.Alt[0] ){
                                allele = 1;
                            }
                            base_q = bam_get_qual(&aln)[querypos];
                        }
                        // the read deletion contain VCF's deletion
                        else if( refAlleleLen != 1 && altAlleleLen == 1 ){

                            if( refAlleleLen != 1 && altAlleleLen == 1){
                                //std::string delSeq = ref_string.substr(ref_pos - 1, align.ol[i] + 1);
                                //std::string refSeq = (*currentVariantIter).second.Ref;
                                //std::string altSeq = (*currentVariantIter).second.Alt;

                                allele = 1;
                                // using this quality to identify indel
                                base_q = -4;
                            }
                            else if ( allele == -1 ) {
                                allele = 0;
                                // using this quality to identify indel
                                base_q = -4;
                            }
                            
                        }
                        
                        if(allele != -1){
                            (*tmpReadResult).variantVec.emplace_back((*currentVariantIter).first, allele, base_q, VariantType::SNP);
                            (*chrMethMap)[variantPos].variantType = VariantType::SNP;
                            currentVariantIter++;
                        }

                    }
                }
            }
            refpos += length;
            ref_pos += length;
        }
        // 3: skipped region from the reference
        else if(cigar_op == 3){
            refpos += length;
            ref_pos += length;
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
    
    if( (*tmpReadResult).variantVec.size() > 0 ){
        (*tmpReadResult).sort();
        readVariantVec.emplace_back((*tmpReadResult));
    }
    delete tmpReadResult;
}

void MethBamParser::exportResult(std::string chrName, std::string chrSquence, int chrLen, std::vector<int> &passPosition, std::ostringstream &modCallResult){

    if(params->outputAllMod){
        for(const auto& pair : *chrMethMap){
            int currPos = pair.first;
            const auto& info = pair.second;
            
            auto posinfoIter = chrMethMap->find(currPos);
            if (posinfoIter == chrMethMap->end()) return;
            std::string infostr = "";
            std::string eachpos;
            std::string samplestr;
            std::string strandinfo;
            std::string ref;
            
            if (chrLen < currPos) return;

            ref = chrSquence.substr(currPos, 1);
            if (ref != "A" && ref != "T" && ref != "C" && ref != "G" &&
                ref != "a" && ref != "t" && ref != "c" && ref != "g") return;

            if (info.strand == 1) strandinfo = "RS=N;";
            else if (info.strand == 0) strandinfo = "RS=P;";
            else return;

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
            samplestr = (*posinfoIter).second.heterstatus + ":" + std::to_string((*posinfoIter).second.methreadcnt) + ":" + std::to_string((*posinfoIter).second.canonreadcnt) + ":" + std::to_string((*posinfoIter).second.depth);

            eachpos = chrName + "\t" + std::to_string((*posinfoIter).first + 1) + "\t" + "." + "\t" + ref + "\t" + "N" + "\t" + "." + "\t" + "PASS" + "\t" + strandinfo + infostr + "\t" + "GT:MD:UD:DP" + "\t" + samplestr + "\n";

            modCallResult<<eachpos;
        }
    }
    else{
        // used to track processed positions, avoid duplicate output
        std::set<int> processedPositions;
        // passPosition is already sorted, no need to sort again
        for(const auto& pos : passPosition){
            // if the position is already processed, skip
            if(processedPositions.count(pos) > 0) {
                continue;
            }

            // process position pos
            auto posinfoIter = chrMethMap->find(pos);
            if (posinfoIter != chrMethMap->end()) {
                std::string infostr = "";
                std::string eachpos;
                std::string samplestr;
                std::string strandinfo;
                std::string ref;
                
                if (chrLen < pos) continue;

                // avoid abnormal REF allele
                ref = chrSquence.substr(pos, 1);
                if (ref != "A" && ref != "T" && ref != "C" && ref != "G" &&
                    ref != "a" && ref != "t" && ref != "c" && ref != "g") continue;

                if (posinfoIter->second.strand == 1) strandinfo = "RS=N;";
                else if (posinfoIter->second.strand == 0) strandinfo = "RS=P;";
                else continue;

                // add modified reads
                if(posinfoIter->second.modReadVec.size() > 0) {
                    infostr += "MR=";
                    for(auto& readName : posinfoIter->second.modReadVec) {
                        infostr += readName + ",";
                    }
                    infostr.back() = ';';
                }
                
                // add non-modified reads
                if(posinfoIter->second.nonModReadVec.size() > 0) {
                    infostr += "NR=";
                    for(auto& readName : posinfoIter->second.nonModReadVec) {
                        infostr += readName + ",";
                    }
                    infostr.back() = ';';
                }
                
                // output all modified positions or only heterozygous positions
                if(params->outputAllMod || posinfoIter->second.heterstatus == "0/1") {
                    int nonmethcnt = posinfoIter->second.canonreadcnt;
                    samplestr = posinfoIter->second.heterstatus + ":" + std::to_string(posinfoIter->second.methreadcnt) + ":" + std::to_string(nonmethcnt) + ":" + std::to_string(posinfoIter->second.depth);
                    
                    eachpos = chrName + "\t" + std::to_string(pos + 1) + "\t" + "." + "\t" + ref + "\t" + "N" + "\t" + "." + "\t" + "PASS" + "\t" + strandinfo + infostr + "\t" + "GT:MD:UD:DP" + "\t" + samplestr + "\n";
                    
                    modCallResult << eachpos;
                }
            }
            processedPositions.insert(pos);
            
            // process position pos+1
            int nextPos = pos + 1;
            auto nextPosinfoIter = chrMethMap->find(nextPos);
            if (nextPosinfoIter != chrMethMap->end() && processedPositions.count(nextPos) == 0) {
                std::string infostr = "";
                std::string eachpos;
                std::string samplestr;
                std::string strandinfo;
                std::string ref;
                
                if (chrLen < nextPos) continue;

                // avoid abnormal REF allele
                ref = chrSquence.substr(nextPos, 1);
                if (ref != "A" && ref != "T" && ref != "C" && ref != "G" &&
                    ref != "a" && ref != "t" && ref != "c" && ref != "g") continue;

                if (nextPosinfoIter->second.strand == 1) strandinfo = "RS=N;";
                else if (nextPosinfoIter->second.strand == 0) strandinfo = "RS=P;";
                else continue;

                // add modified reads
                if(nextPosinfoIter->second.modReadVec.size() > 0) {
                    infostr += "MR=";
                    for(auto& readName : nextPosinfoIter->second.modReadVec) {
                        infostr += readName + ",";
                    }
                    infostr.back() = ';';
                }
                
                // add non-modified reads
                if(nextPosinfoIter->second.nonModReadVec.size() > 0) {
                    infostr += "NR=";
                    for(auto& readName : nextPosinfoIter->second.nonModReadVec) {
                        infostr += readName + ",";
                    }
                    infostr.back() = ';';
                }
                
                // output all modified positions or only heterozygous positions
                if(params->outputAllMod || nextPosinfoIter->second.heterstatus == "0/1") {
                    int nonmethcnt = nextPosinfoIter->second.canonreadcnt;
                    samplestr = nextPosinfoIter->second.heterstatus + ":" + std::to_string(nextPosinfoIter->second.methreadcnt) + ":" + std::to_string(nonmethcnt) + ":" + std::to_string(nextPosinfoIter->second.depth);
                    
                    eachpos = chrName + "\t" + std::to_string(nextPos + 1) + "\t" + "." + "\t" + ref + "\t" + "N" + "\t" + "." + "\t" + "PASS" + "\t" + strandinfo + infostr + "\t" + "GT:MD:UD:DP" + "\t" + samplestr + "\n";
                    
                    modCallResult << eachpos;
                }
                processedPositions.insert(nextPos);
            }
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

void MethBamParser::judgeMethGenotype(std::string chrName, std::vector<ReadVariant> &readVariantVec, std::vector<ReadVariant> &modReadVariantVec){
    //Find all methylation sites and determine their genotypes
    for(std::map<int, MethPosInfo>::iterator chrmethmapIter = chrMethMap->begin(); chrmethmapIter != chrMethMap->end(); chrmethmapIter++){

        // Determine methylation genotype
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
        
        // Determine methylation genotype for each position individually
        if(heterRatio >= params->heterRatio && noiseRatio <= params->noiseRatio){
            (*chrmethmapIter).second.heterstatus = "0/1";
        }
        else if(methcnt >= nonmethcnt){
            (*chrmethmapIter).second.heterstatus = "1/1";
        }
        else{
            (*chrmethmapIter).second.heterstatus = "0/0";
        }
    }
    
    std::set<int> positionPairs; 

    // Only process cases with forward-reverse strand pairs
    for(auto iter = chrMethMap->begin(); iter != chrMethMap->end(); iter++){        
        if(iter->second.strand == 0 && iter->second.variantType == VariantType::MOD){
            int currPos = iter->first;
            int nextPos = currPos + 1; 
            
            auto nextIter = chrMethMap->find(nextPos);
            // Check that the next position exists, is on the reverse strand, and is a methylation site
            if(nextIter != chrMethMap->end() && nextIter->second.strand == 1 && nextIter->second.variantType == VariantType::MOD){
                // Calculate methylation statistics after combining
                float totalMethCnt = iter->second.methreadcnt + nextIter->second.methreadcnt;
                float totalNonMethCnt = iter->second.canonreadcnt + nextIter->second.canonreadcnt;
                float totalDepth = iter->second.depth + nextIter->second.depth;
                float totalNoiseCnt = totalDepth - totalMethCnt - totalNonMethCnt;
                
                if(std::max(totalMethCnt, totalNonMethCnt) == 0){
                    continue;
                }
                
                float combinedHeterRatio = std::min(totalMethCnt, totalNonMethCnt) / std::max(totalMethCnt, totalNonMethCnt);
                float combinedNoiseRatio = totalNoiseCnt / totalDepth;
                
                // Determine methylation genotype based on combined data
                std::string combinedStatus;
                if(combinedHeterRatio >= params->heterRatio && combinedNoiseRatio <= params->noiseRatio){
                    combinedStatus = "0/1";
                    positionPairs.insert(currPos);
                }
                else if(totalMethCnt >= totalNonMethCnt){
                    combinedStatus = "1/1";
                }
                else{
                    combinedStatus = "0/0";
                }                
                // Update the status of forward and reverse strand positions
                iter->second.heterstatus = combinedStatus;
                nextIter->second.heterstatus = combinedStatus;
            }
        }
    }

    // Handle all variants uniformly to avoid duplication
    for(auto& read : readVariantVec){
        ReadVariant newRead;
        newRead.read_name = read.read_name;
        newRead.is_reverse = read.is_reverse;
        
        for(auto& variant : read.variantVec){
            int pos = variant.position;
            // Check if the variant is methylation-related and marked as heterogeneous in chrMethMap
            if(variant.type == VariantType::MOD) {
                if(chrMethMap->find(pos)->second.strand == 0){
                    auto pairIter = positionPairs.find(pos);
                    if(pairIter != positionPairs.end()){
                        newRead.variantVec.emplace_back(pos, variant.allele, variant.quality, VariantType::MOD);
                    }
                }
                else if(chrMethMap->find(pos)->second.strand == 1){
                    auto pairIter = positionPairs.find(pos-1);
                    if(pairIter != positionPairs.end()){
                        newRead.variantVec.emplace_back(pos-1, variant.allele, variant.quality, VariantType::MOD);
                    }
                }
            } 
            else if(variant.type == VariantType::SNP) {
                newRead.variantVec.emplace_back(variant);
            }
        }
        if(!newRead.variantVec.empty()){
            modReadVariantVec.push_back(newRead);
        }
        newRead.variantVec.clear();
        newRead.variantVec.shrink_to_fit();
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
    nodeInfo = new std::map<int, std::map<std::string, VariantType>>;
    edgeList = new std::map<int,VariantEdge*>;
    forwardModNode = new std::map<int,ReadBaseMap*>;
    reverseModNode = new std::map<int,ReadBaseMap*>;
}

MethylationGraph::~MethylationGraph(){
}

void MethylationGraph::addEdge(std::vector<ReadVariant> &in_readVariant, std::string chrName){
    readVariant = &in_readVariant;
    int readCount = 0;
    int vectorSize = 0;
    // iter all read
    for(std::vector<ReadVariant>::iterator readIter = in_readVariant.begin() ; readIter != in_readVariant.end() ; readIter++ ){
        ReadVariant tmpRead;
        vectorSize += (*readIter).variantVec.size();
        readCount++;
        
        for( auto variant : (*readIter).variantVec ){           
            auto nodeIter = nodeInfo->find(variant.position);
            
            if( nodeIter == nodeInfo->end() ){
                (*nodeInfo)[variant.position] = std::map<std::string, VariantType>();
            }
            (*nodeInfo)[variant.position][(*readIter).read_name] = variant.type;
        }

        for (auto variant1Iter = (*readIter).variantVec.begin(); variant1Iter != (*readIter).variantVec.end(); ++variant1Iter) {
            auto variant2Iter = std::next(variant1Iter);
            int searchCount = 0;
            while (variant2Iter != (*readIter).variantVec.end() && searchCount < 50) {
                if (!(variant1Iter->type == VariantType::SNP && variant2Iter->type == VariantType::SNP) ) {

                    int pos = variant1Iter->position;
                    auto edgeIter = edgeList->find(pos);
                    if (edgeIter == edgeList->end()) {
                        (*edgeList)[pos] = new VariantEdge(pos);
                    }

                    if (variant1Iter->allele == 0) {
                        (*edgeList)[pos]->ref->addSubEdge(variant1Iter->quality, *variant2Iter, readIter->read_name, 0, 1);
                    } 
                    else if (variant1Iter->allele == 1) {
                        (*edgeList)[pos]->alt->addSubEdge(variant1Iter->quality, *variant2Iter, readIter->read_name, 0, 1);
                    }
                }
                ++searchCount;
                ++variant2Iter;
            }
        }
    }
}

void MethylationGraph::connectResults(std::string chrName, std::vector<int> &passPosition, bool hasValidSnpData){
    std::set<int> strongMethylationPoints;
    std::set<int> weakMethylationPoints;
    std::set<int> weakMethylationPoints2;
    std::set<int> addedPositions;
    std::set<int> addedPositions2;
    std::vector<int> prepassPosition;
    std::vector<int> hasConnect;

    // If no valid SNP data, skip the first pass and directly add all methylation sites to strongMethylationPoints
    if (!hasValidSnpData) {
        for (auto nodeIter = nodeInfo->begin(); nodeIter != nodeInfo->end(); ++nodeIter) {
            int currPos = nodeIter->first;
            if(checkVariantType(currPos) == VariantType::MOD){
                strongMethylationPoints.insert(currPos);
            }
        }
    }
    else {
        // First pass: Identify methylation points with strong SNP connections
        for (auto nodeIter = nodeInfo->begin(); nodeIter != nodeInfo->end(); ++nodeIter) {
            int currPos = nodeIter->first;
            auto nextNodeIter = std::next(nodeIter, 1);
            int searchCount = 0;
            if( nextNodeIter == nodeInfo->end() ){
                break;
            }

            auto edgeIter = edgeList->find(currPos);
            if (edgeIter == edgeList->end()) 
                continue;

            if(checkVariantType(currPos) == 0){
                auto searchNodeIter = nextNodeIter;
                // Continue searching until we find a SNP or reach the end of nodeInfo
                while (searchNodeIter != nodeInfo->end() && searchCount < params->connectAdjacent) {
                    std::pair<int, int> tmp = edgeIter->second->findNumberOfRead(searchNodeIter->first);
                    int totalConnectReads = tmp.first + tmp.second;
                    float minimumConnection = std::max(((*nodeIter).second.size() + (*searchNodeIter).second.size())/ 4.0, 6.0);
                    if(totalConnectReads <= minimumConnection){
                        break;
                    }
                    if(checkVariantType(searchNodeIter->first) == VariantType::SNP) {     
                        double majorRatio = (double)std::max(tmp.first, tmp.second) / (double)(tmp.first + tmp.second);
                        hasConnect.push_back(currPos);
                        if (majorRatio >= params->connectConfidence && totalConnectReads > minimumConnection && strongMethylationPoints.count(currPos) == 0) {
                            strongMethylationPoints.insert(currPos);
                            break;
                        }
                    }
                    ++searchNodeIter;
                    ++searchCount;
                }
                if(std::find(hasConnect.begin(), hasConnect.end(), currPos) == hasConnect.end()){
                    weakMethylationPoints.insert(currPos);
                }
            }
            else if(checkVariantType(currPos) == VariantType::SNP){
                auto searchNodeIter = nextNodeIter;
                prepassPosition.push_back(currPos);
                while(searchNodeIter != nodeInfo->end()) {
                    std::pair<int, int> tmp = edgeIter->second->findNumberOfRead(searchNodeIter->first);
                    int totalConnectReads = tmp.first + tmp.second;
                    float minimumConnection = std::max(((*nodeIter).second.size() + (*searchNodeIter).second.size())/ 4.0, 6.0);
                    if(totalConnectReads <= minimumConnection){
                        break;
                    }
                    if (checkVariantType(searchNodeIter->first) == VariantType::MOD) {
                        double majorRatio = (double)std::max(tmp.first, tmp.second) / (double)(tmp.first + tmp.second);
                        hasConnect.push_back(searchNodeIter->first);
                        if (majorRatio >= params->connectConfidence && totalConnectReads > minimumConnection && strongMethylationPoints.count(nextNodeIter->first) == 0) {
                            strongMethylationPoints.insert(nextNodeIter->first);
                        } 
                    }
                    ++searchNodeIter;      
                    ++searchCount;
                }
            }
        }
    } 
    hasConnect.clear();

    // Second pass: Evaluate connections between strong methylation points
    for(auto it1 = strongMethylationPoints.begin(); it1 != strongMethylationPoints.end(); ++it1){
        int pos1 = *it1;
        auto it2 = std::next(it1, 1);
        auto searchNodeIter = it2;
        int searchCount = 0;
        auto edgeIter = edgeList->find(pos1);
        if (edgeIter == edgeList->end()) 
            continue;

        while(searchNodeIter != strongMethylationPoints.end() && searchCount < params->connectAdjacent) {
            int pos2 = *searchNodeIter;
            std::pair<int, int> tmp = edgeIter->second->findNumberOfRead(pos2);
            int totalConnectReads = tmp.first + tmp.second;
            float minimumConnection = std::max(((*nodeInfo)[pos1].size() + (*nodeInfo)[pos2].size())/ 4.0, 6.0);
            double majorRatio = (double)std::max(tmp.first, tmp.second) / (double)(tmp.first + tmp.second);

            if(totalConnectReads <= minimumConnection){
                break;
            }

            if (majorRatio >= params->connectConfidence && totalConnectReads > minimumConnection) {
                if(addedPositions.count(pos1) == 0) {
                    prepassPosition.push_back(pos1);
                    addedPositions.insert(pos1);
                    // Only add positions to weakMethylationPoints if there is valid SNP data (for third pass)
                    if(hasValidSnpData) {
                        weakMethylationPoints.insert(pos1);
                    }
                }
                if(addedPositions.count(pos2) == 0) {
                    prepassPosition.push_back(pos2);
                    addedPositions.insert(pos2);
                    // Only add positions to weakMethylationPoints if there is valid SNP data (for third pass)
                    if(hasValidSnpData) {
                        weakMethylationPoints.insert(pos2);
                    }
                }
            }
            ++searchNodeIter;
            ++searchCount;
        }
    }
    strongMethylationPoints.clear();

    //third pass: evaluate connections between weak methylation points
    for(int i = 0; i < params->iterCount; i++){
        if (hasValidSnpData) {
            // Use alternating sets for each iteration
            auto& currentWeakPoints = (i % 2 == 0) ? weakMethylationPoints : weakMethylationPoints2;
            auto& nextWeakPoints = (i % 2 == 0) ? weakMethylationPoints2 : weakMethylationPoints;
            auto& currentAddedPositions = (i % 2 == 0) ? addedPositions : addedPositions2;
            auto& nextAddedPositions = (i % 2 == 0) ? addedPositions2 : addedPositions;
            
            // Clear the target set for this iteration
            nextWeakPoints.clear();
            nextAddedPositions.clear();
            
            for(auto it1 = currentWeakPoints.begin(); it1 != currentWeakPoints.end(); ++it1){
                int currPos = *it1;
                auto nextIter = it1;
                int nextSearchCount = 0;
                bool isAdded = false;
                auto edgeIter = edgeList->find(currPos);
                if (edgeIter == edgeList->end()) 
                    continue;
                while(++nextIter != currentWeakPoints.end() && nextSearchCount < params->connectAdjacent) {
                    int nextPos = *nextIter;
                    if(currentAddedPositions.count(nextPos) == 0 && currentAddedPositions.count(currPos) == 0){
                        ++nextSearchCount;
                        continue;
                    }
                    isAdded = true;
                    std::pair<int, int> tmp = edgeIter->second->findNumberOfRead(nextPos);
                    int totalConnectReads = tmp.first + tmp.second;
                    float minimumConnection = std::max(((*nodeInfo)[currPos].size() + (*nodeInfo)[nextPos].size())/ 4.0, 6.0);
                    double majorRatio = (double)std::max(tmp.first, tmp.second) / (double)(tmp.first + tmp.second);
                    if(totalConnectReads <= minimumConnection){
                        break;
                    }
                    if (majorRatio >= params->connectConfidence && totalConnectReads > minimumConnection ) {
                        if(std::find(prepassPosition.begin(), prepassPosition.end(), currPos) == prepassPosition.end()){
                            prepassPosition.push_back(currPos);
                            nextWeakPoints.insert(currPos);
                            nextAddedPositions.insert(currPos);
                        }
                        if(std::find(prepassPosition.begin(), prepassPosition.end(), nextPos) == prepassPosition.end()){
                            prepassPosition.push_back(nextPos);
                            nextWeakPoints.insert(nextPos);
                            nextAddedPositions.insert(nextPos);
                        }
                    }
                    ++nextSearchCount;
                }
                if(!isAdded){
                    nextWeakPoints.insert(currPos);
                }
            }
        }
    }
    weakMethylationPoints.clear();
    weakMethylationPoints2.clear();
    addedPositions.clear();
    addedPositions2.clear();

    // Ensure passPosition is sorted by position
    std::sort(prepassPosition.begin(), prepassPosition.end());
    // Fourth step: Filter positions that do not have good connections to both neighbors
    for (size_t i = 0; i < prepassPosition.size(); ++i) {
        int pos = prepassPosition[i];
        bool hasGoodPrevConnection = false;
        bool hasGoodNextConnection = false;
        if(nodeInfo->find(pos) != nodeInfo->end()){
            if(checkVariantType(pos) == VariantType::SNP){
                continue;
            }
        }
                
        // Check connection with previous position
        if (i > 0) {
            int prevPos = prepassPosition[i-1];
            auto edgeIter = edgeList->find(prevPos);
            if (edgeIter == edgeList->end()) {
                hasGoodPrevConnection = true;
                continue;
            }
            std::pair<int, int> tmp = edgeIter->second->findNumberOfRead(pos);
            int totalConnectReads = tmp.first + tmp.second;
            double majorRatio = (double)std::max(tmp.first, tmp.second) / (double)(tmp.first + tmp.second);
            if(totalConnectReads != 0){
                
                if (majorRatio >= params->connectConfidence && totalConnectReads >= 6 ) {
                    hasGoodPrevConnection = true;
                }
            }
        }
        
        // Check connection with next position
        if (i < prepassPosition.size() - 1 && hasGoodPrevConnection) {
            int nextPos = prepassPosition[i+1];
            auto edgeIter = edgeList->find(pos);
            if (edgeIter == edgeList->end()) {
                hasGoodNextConnection = true;
                continue;
            }
            std::pair<int, int> tmp = edgeIter->second->findNumberOfRead(nextPos);
            int totalConnectReads = tmp.first + tmp.second;
            double majorRatio = (double)std::max(tmp.first, tmp.second) / (double)(tmp.first + tmp.second);
            if(totalConnectReads != 0){
                if (majorRatio >= params->connectConfidence && totalConnectReads >= 6 ) {
                    hasGoodNextConnection = true;
                }
            }
        }
        
        // Only keep positions with good connections to both neighbors (or edge positions)
        if ( hasGoodNextConnection || i == 0 || i == prepassPosition.size() - 1) {
            passPosition.push_back(pos);
        }
    }
    prepassPosition.clear();
}

int MethylationGraph::checkVariantType(int position){
    auto nodeIter = nodeInfo->find(position);
    if (nodeIter != nodeInfo->end()) {
        for (const auto& nodeType : nodeIter->second) {
            if(nodeType.second == VariantType::MOD){
                return VariantType::MOD; // Return true if any type is MOD
            }
            else if(nodeType.second == VariantType::SNP){
                return VariantType::SNP; // Return true if any type is SNP
            }
            else if(nodeType.second == VariantType::INDEL){
                return VariantType::INDEL; // Return true if any type is INDEL
            }
            else if(nodeType.second == VariantType::SV){
                return VariantType::SV; // Return true if any type is SV
            }
            else{
                return -1; // Return -1 if no variant type is found
            }
        }
    }
    return -1; // Return -1 if the position is not in the nodeInfo
}

void MethylationGraph::destroy(){

    for( auto edgeIter = edgeList->begin() ; edgeIter != edgeList->end() ; edgeIter++ ){
        edgeIter->second->ref->destroy();
        edgeIter->second->alt->destroy();
        delete edgeIter->second->ref;
        delete edgeIter->second->alt;
    }
    
    for( auto nodeIter = forwardModNode->begin() ; nodeIter != forwardModNode->end() ; nodeIter++ ){
        delete nodeIter->second;
    }
    
    for( auto nodeIter = reverseModNode->begin() ; nodeIter != reverseModNode->end() ; nodeIter++ ){
        delete nodeIter->second;
    }
    
    delete nodeInfo;
    delete edgeList;
}

MethSnpParser::MethSnpParser(ModCallParameters &in_params):commandLine(false), hasSnpData(false){
    chrVariant = new std::map<std::string, std::map<int, RefAlt> >;
    
    params = &in_params;
    
    // Check if SNP file is provided
    if(params->snpFile.empty()) {
        std::cerr << "No SNP file provided, running without SNP variants.";
        return; 
    }
    
    // open vcf file
    htsFile * inf = bcf_open(params->snpFile.c_str(), "r");
    if(inf == nullptr) {
        std::cerr << "Warning: Could not open SNP file " << params->snpFile << ", running without SNP variants.";
        return;
    }
    
    // read header
    bcf_hdr_t *hdr = bcf_hdr_read(inf);
    if(hdr == nullptr) {
        std::cerr << "Warning: Could not read SNP file header, running without SNP variants.";
        bcf_close(inf);
        return;
    }
    // counters
    int nseq = 0;
    // report names of all the sequences in the VCF file
    const char **seqnames = NULL;
    // chromosome idx and name
    seqnames = bcf_hdr_seqnames(hdr, &nseq);
    // store chromosome
    for (int i = 0; i < nseq; i++) {
        // bcf_hdr_id2name is another way to get the name of a sequence
        chrName.emplace_back(seqnames[i]);
    }
    // set all sample string
    std::string allSmples = "-";
    // limit the VCF data to the sample name passed in
    int is_file = bcf_hdr_set_samples(hdr, allSmples.c_str(), 0);
    if( is_file != 0 ){
        std::cout << "error or a positive integer if the list contains samples not present in the VCF header\n";
    }

    // struct for storing each record
    bcf1_t *rec = bcf_init();
    int ngt_arr = 0;
    int ngt = 0;
    int *gt     = NULL;
    
    // loop vcf line 
    while (bcf_read(inf, hdr, rec) == 0) {
        // snp
        if (bcf_is_snp(rec)) {
            ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);

            if(ngt<0){
                std::cerr<< "pos " << rec->pos << " missing GT value" << "\n";
                exit(1);
            }
            
            // just phase hetero SNP
            if ( (gt[0] == 2 && gt[1] == 4) || // 0/1
                 (gt[0] == 4 && gt[1] == 2) || // 1/0
                 (gt[0] == 2 && gt[1] == 5) || // 0|1
                 (gt[0] == 4 && gt[1] == 3)    // 1|0 
                ) {
                
                // get chromosome string
                std::string chr = seqnames[rec->rid];
                // position is 0-base
                int variantPos = rec->pos;
                // get r alleles
                RefAlt tmp;
                tmp.Ref = rec->d.allele[0]; 
                tmp.Alt = rec->d.allele[1];
                
                //prevent the MAVs calling error which makes the GT=0/1
                if ( rec->d.allele[1][2] != '\0' ){
                    continue;
                }
                
                // record 
                (*chrVariant)[chr][variantPos] = tmp;
            }
        }
    }
    
    // Clean up resources
    if(gt) free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(inf);
    
    // Set successful loading flag
    hasSnpData = true;
    std::cerr << "Successfully loaded SNP variants from " << params->snpFile << "\n";
}

void MethSnpParser::parserProcess(std::string &input) {
}

MethSnpParser::~MethSnpParser() {
    delete chrVariant;
}

bool MethSnpParser::hasValidSnpData() const {
    return hasSnpData;
}


void MethSnpParser::writeLine(std::string &input, bool &ps_def, std::ofstream &resultVcf, PhasingResult &phasingResult) {
    //header
    if( input.substr(0, 2) == "##" ){
        // avoid double definition
        if( input.substr(0, 16) == "##FORMAT=<ID=PS," ){
            ps_def = true;
        }
        resultVcf << input << "\n";
    }
    else if ( input.substr(0, 6) == "#CHROM" || input.substr(0, 6) == "#chrom" ){
        // format line 
        if( commandLine == false ){
            if(ps_def == false){
                resultVcf <<  "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">\n";
                ps_def = true;
            }
            resultVcf << "##longphaseVersion=" << params->version << "\n";
            resultVcf << "##commandline=\"" << params->command << "\"\n";
            commandLine = true;
        }
        resultVcf << input << "\n";
    }
    else{
        std::istringstream iss(input);
        std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

        if( fields.size() == 0 )
            return;

        int pos = std::stoi( fields[1] );
        int posIdx = pos - 1 ;
        std::string key = fields[0] + "_" + std::to_string(posIdx);

        PhasingResult::iterator psElementIter =  phasingResult.find(key);

        // PS flag already exist, erase PS info
        if( fields[8].find("PS")!= std::string::npos ){

            // find PS flag
            int colon_pos = 0;
            int ps_pos = fields[8].find("PS");
            for(int i =0 ; i< ps_pos ; i++){
                if(fields[8][i]==':')
                    colon_pos++;
            }
                                
            // erase PS flag
            if( fields[8].find(":",ps_pos+1) != std::string::npos ){
                fields[8].erase(ps_pos, 3);
            }
            else{
                fields[8].erase(ps_pos - 1, 3);
            }
                            
            // find PS value start
            int current_colon = 0;
            int ps_start = 0;
            for(unsigned int i =0; i < fields[9].length() ; i++){
                if( current_colon >= colon_pos )
                    break;
                if(fields[9][i]==':')
                    current_colon++;  
                ps_start++;
            }
                            
            // erase PS value
            if( fields[9].find(":",ps_start+1) != std::string::npos ){
                int ps_end_pos = fields[9].find(":",ps_start+1);
                fields[9].erase(ps_start, ps_end_pos - ps_start + 1);
            }
            else{
                fields[9].erase(ps_start-1, fields[9].length() - ps_start + 1);
            }
        }
        // reset GT flag
        if( fields[8].find("GT")!= std::string::npos ){
            // find GT flag
            int colon_pos = 0;
            int gt_pos = fields[8].find("GT");
            for(int i =0 ; i< gt_pos ; i++){
                if(fields[8][i]==':')
                    colon_pos++;
            }
            // find GT value start
            int current_colon = 0;
            int modify_start = 0;
            for(unsigned int i =0; i < fields[9].length() ; i++){
                if( current_colon >= colon_pos )
                    break;
                if(fields[9][i]==':')
                    current_colon++;  
                modify_start++;
            }
            if(fields[9][modify_start+1] == '|'){
                // direct modify GT value
                if(fields[9][modify_start] > fields[9][modify_start+2]){
                    fields[9][modify_start+1] = fields[9][modify_start];
                    fields[9][modify_start] = fields[9][modify_start+2];
                    fields[9][modify_start+2] =fields[9][modify_start+1];
                }
                fields[9][modify_start+1] = '/';
            }   
        }

        // Check if the variant is extracted from this VCF
        auto posIter = (*chrVariant)[fields[0]].find(posIdx);
        
        // this pos is phase
        if( psElementIter != phasingResult.end() && posIter != (*chrVariant)[fields[0]].end() ){
            // add PS flag and value
            fields[8] = fields[8] + ":PS";
            fields[9] = fields[9] + ":" + std::to_string((*psElementIter).second.block);
                        
            // find GT flag
            int colon_pos = 0;
            int gt_pos = fields[8].find("GT");
            for(int i =0 ; i< gt_pos ; i++){
                if(fields[8][i]==':')
                    colon_pos++;
                }
            // find GT value start
            int current_colon = 0;
            int modify_start = 0;
            for(unsigned int i =0; i < fields[9].length() ; i++){
                if( current_colon >= colon_pos )
                    break;
                if(fields[9][i]==':')
                    current_colon++;  
                modify_start++;
            }
            // direct modify GT value
            fields[9][modify_start] = (*psElementIter).second.RAstatus[0];
            fields[9][modify_start+1] = '|';
            fields[9][modify_start+2] = (*psElementIter).second.RAstatus[2];
        }
        // this pos has not been phased
        else{
            // add PS flag and value
            fields[8] = fields[8] + ":PS";
            fields[9] = fields[9] + ":.";
        }
                    
        for(std::vector<std::string>::iterator fieldIter = fields.begin(); fieldIter != fields.end(); ++fieldIter){
            if( fieldIter != fields.begin() )
                resultVcf<< "\t";
            resultVcf << (*fieldIter);
        }
        resultVcf << "\n";
    }
}

std::map<int, RefAlt> MethSnpParser::getVariants(std::string chrName) {
    std::map<int, RefAlt> targetVariants;
    std::map<std::string, std::map<int, RefAlt> >::iterator chrIter = chrVariant->find(chrName);
    
    if( chrIter != chrVariant->end() )
        targetVariants = (*chrIter).second;
    
    return targetVariants;
}

std::map<int, RefAlt> MethSnpParser::getVariants_markindel(std::string chrName, const std::string &ref) {
    std::map<int, RefAlt> targetVariants;
    std::map<std::string, std::map<int, RefAlt> >::iterator chrIter = chrVariant->find(chrName);

    //Mark the indel which lies in the tandem repeat
    if (chrIter != chrVariant->end()) {
        // Accessing the inner map using chrIter->second and iterating over it
        for ( auto innerIter = chrIter->second.begin(); innerIter != chrIter->second.end(); innerIter++ ) {
            int variant_pos = innerIter->first;      // Accessing the variant_pos of the inner map
            RefAlt variant_info = innerIter->second; // Accessing the variant_info of the inner map
            int ref_pos = variant_pos ;
	        bool danger = false ;

            std::string repeat = ref.substr(ref_pos + 1, 2) ; // get the 2 words string in the reference which behind the indel position
            int i = 0 ;
	        //check if there has 2 words tandem repeat in the reference
            while ( i < 5 && (variant_info.Ref.length()>1 ||variant_info.Alt.length()>1) /*&& repeat[0]!=repeat[1]*/ ) {
                if ( repeat[0] != ref[ref_pos+1] || repeat[1] != ref[ref_pos+2] ) {
                    break ;
                }

                ref_pos = ref_pos + 2 ;
                i++ ;
            }

	    // set the dnager to true if the repeat word repeats at least five times
	    if ( i == 5 ) danger = true ;

            innerIter->second.is_danger = danger ;

        }

        targetVariants = (*chrIter).second;
    }
    else {
        // Handle the case where chrName doesn't exist in chrVariant
    }

    return targetVariants;
}

std::vector<std::string> MethSnpParser::getChrVec() {
    return chrName;
}

int MethSnpParser::getLastSNP(std::string chrName) {
    std::map<std::string, std::map< int, RefAlt > >::iterator chrVariantIter = chrVariant->find(chrName);
    // this chromosome not exist in this file. 
    if(chrVariantIter == chrVariant->end())
        return -1;
    // get last SNP
    std::map< int, RefAlt >::reverse_iterator lastVariantIter = (*chrVariantIter).second.rbegin();
    // there are no SNPs on this chromosome
    if(lastVariantIter == (*chrVariantIter).second.rend())
        return -1;
    return (*lastVariantIter).first;
}

void MethSnpParser::writeResult(PhasingResult phasingResult) {
    if( params->snpFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressInput(params->snpFile, params->resultPrefix+".vcf", phasingResult);
    }
    else if( params->snpFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressInput(params->snpFile, params->resultPrefix+".vcf", phasingResult);
    }
    return;
}
