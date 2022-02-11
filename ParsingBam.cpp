#include "ParsingBam.h"

#include <string.h>
#include <sstream>


// vcf parser modify from 
// http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html
// bam parser modify from 
// https://github.com/whatshap/whatshap/blob/9882248c722b1020321fea6e9491c1cc5b75354b/whatshap/variants.py

FastaParser::FastaParser(std::string fastaFile ,  std::vector<std::string> chrName , std::vector<int> last_pos):
    fastaFile(fastaFile),
    chrName(chrName),
    last_pos(last_pos)
{
    if(fastaFile==""){
        for(std::vector<std::string>::iterator iter = chrName.begin() ; iter != chrName.end() ; iter++)
            chrString.insert(std::make_pair( (*iter) , ""));
        return;
    }
    
    faidx_t *fai = NULL;
    fai = fai_load(fastaFile.c_str());
    // iterating all chr
    for(std::vector<std::string>::iterator iter = chrName.begin() ; iter != chrName.end() ; iter++){
        int index = iter - chrName.begin();

        // ref_len is a return value that is length of retrun string
        int ref_len = 0;
        // read file
        std::string chr_info(faidx_fetch_seq(fai , (*iter).c_str() , 0 ,last_pos.at(index)+5 , &ref_len));
        if(ref_len == 0){
            std::cout<<"nothing in reference file \n";
        }
        // insert to map
        chrString.insert(std::make_pair( (*iter) , chr_info));

    }
}
FastaParser::~FastaParser(){
 
}

VariantParser::VariantParser(std::string variantFile):variantFile(variantFile){
    // open vcf file
    htsFile * inf = bcf_open(variantFile.c_str(), "r");
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
    if( is_file != 0 ){
        std::cout << "error or a positive integer if the list contains samples not present in the VCF header\n";
    }
    
    //if (bcf_hdr_nsamples(hdr) != 1) {
    //    fprintf(stderr, "ERROR: please limit to a single sample\n");
    //}

    // struct for storing each record
    bcf1_t *rec = bcf_init();
    int ngt_arr = 0;
    int ngt = 0;
    int *gt     = NULL;
    
    //int ps_arr = 0;
    //int nps = 0;
    //int *ps    = NULL;
       
    // loop vcf line 
    while (bcf_read(inf, hdr, rec) == 0) {
        // this function directly looking Ref and Alt string length
        if (bcf_is_snp(rec)) {
            ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
            /*
            nps = bcf_get_format_int32(hdr, rec, "PS", &ps, &ps_arr);
            if(nps>0){
                std::cerr<< "pos " << rec->pos << " already phased\n";
                std::cerr<< "please check vcf file.\n";
                exit(1);
            }
            */
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
                // record 
                chrVariant[chr][variantPos] = tmp;
            }
        } 
    }
}
VariantParser::~VariantParser(){
    
}

std::map<int, RefAlt> VariantParser::getVariants(std::string chrName){
    std::map<int, RefAlt> targetVariants;
    std::map<std::string, std::map<int, RefAlt> >::iterator chrIter = chrVariant.find(chrName);
    
    if( chrIter != chrVariant.end() )
        targetVariants = (*chrIter).second;
    
    return targetVariants;
}

std::vector<std::string> VariantParser::getChrVec(){
    return chrName;
}

bool VariantParser::findChromosome(std::string chrName){
    std::map<std::string, std::map< int, RefAlt > >::iterator chrVariantIter = chrVariant.find(chrName);
    // this chromosome not exist in this file. 
    if(chrVariantIter == chrVariant.end())
        return false;
    return true;
}

int VariantParser::getLastSNP(std::string chrName){
    std::map<std::string, std::map< int, RefAlt > >::iterator chrVariantIter = chrVariant.find(chrName);
    // this chromosome not exist in this file. 
    if(chrVariantIter == chrVariant.end())
        return -1;
    // get last SNP
    std::map< int, RefAlt >::reverse_iterator lastVariantIter = (*chrVariantIter).second.rbegin();
    // there are no SNPs on this chromosome
    if(lastVariantIter == (*chrVariantIter).second.rend())
        return -1;
    return (*lastVariantIter).first;
}

void VariantParser::writeResult(std::string resultPrefix, PhasingResult phasingResult){

    if( variantFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressInput(resultPrefix,phasingResult);
    }
    else if( variantFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressInput(resultPrefix,phasingResult);
    }
    return;
}

void VariantParser::compressInput(std::string resultPrefix, PhasingResult phasingResult){
    
    gzFile file = gzopen(variantFile.c_str(), "rb");
    std::ofstream resultVcf(resultPrefix+".vcf");
    
    if(!resultVcf.is_open()){
        std::cout<< "Fail to open write file: " << resultPrefix+".vcf" << "\n";
    }
    else if(!file){
        std::cout<< "Fail to open vcf: " << variantFile << "\n";
    }
    else{
        bool ps_def = false;
              
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
                writeLine(input, ps_def, resultVcf, phasingResult);
            }
            // any trailing data in [eol, end) now is a partial line
            offset = std::copy(cur, end, buffer);
        }
        gzclose (file);
    }
}

void VariantParser::unCompressInput(std::string resultPrefix, PhasingResult phasingResult){
    std::ifstream originVcf(variantFile);
    std::ofstream resultVcf(resultPrefix+".vcf");
    
    if(!resultVcf.is_open()){
        std::cout<< "Fail to open write file: " << resultPrefix+".vcf" << "\n";
    }
    else if(!originVcf.is_open()){
        std::cout<< "Fail to open vcf: " << variantFile << "\n";
    }
    else{
        bool ps_def = false;
        std::string input;
        while(! originVcf.eof() ){
            std::getline(originVcf, input);
            writeLine(input, ps_def, resultVcf, phasingResult);
        }
    }
}
void VariantParser::writeLine(std::string &input, bool &ps_def, std::ofstream &resultVcf, PhasingResult &phasingResult){
    // header
    if( input.find("#")!= std::string::npos){
        // avoid double definition
        if( input.find("ID=PS")!= std::string::npos){
            ps_def = true;
        }
        // format line 
        if(input.find("##")== std::string::npos && ps_def == false){
            resultVcf <<  "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">\n";
        }
        resultVcf << input << "\n";
    }
    else if( input.find("#")== std::string::npos ){
        std::istringstream iss(input);
        std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

        if( fields.size() == 0 )
            return;

        int pos = std::stoi( fields[1] );
        int posIdx = pos - 1 ;
        std::string key = fields[0] + "_" + std::to_string(posIdx);

        PhasingResult::iterator psElementIter =  phasingResult.find(key);
        std::map<int, RefAlt>::iterator posIter = chrVariant[fields[0]].find(posIdx);

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
            int modifu_start = 0;
            for(unsigned int i =0; i < fields[9].length() ; i++){
                if( current_colon >= colon_pos )
                    break;
                if(fields[9][i]==':')
                    current_colon++;  
                modifu_start++;
            }
            if(fields[9][modifu_start+1] == '|'){
                // direct modify GT value
                if(fields[9][modifu_start] > fields[9][modifu_start+2]){
                    fields[9][modifu_start+1] = fields[9][modifu_start];
                    fields[9][modifu_start] = fields[9][modifu_start+2];
                    fields[9][modifu_start+2] =fields[9][modifu_start+1];
                }
                fields[9][modifu_start+1] = '/';
            }   
        }

        // this pos is phase
        if( psElementIter != phasingResult.end() && posIter != chrVariant[fields[0]].end()){
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
            int modifu_start = 0;
            for(unsigned int i =0; i < fields[9].length() ; i++){
                if( current_colon >= colon_pos )
                    break;
                if(fields[9][i]==':')
                    current_colon++;  
                modifu_start++;
            }
            // direct modify GT value
            fields[9][modifu_start] = (*psElementIter).second.RAstatus[0];
            fields[9][modifu_start+1] = '|';
            fields[9][modifu_start+2] = (*psElementIter).second.RAstatus[2];
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


bool VariantParser::findSNP(std::string chr, int posistion){
    std::map<std::string, std::map<int, RefAlt> >::iterator chrIter = chrVariant.find(chr);
    // empty chromosome
    if( chrIter == chrVariant.end() )
        return false;
    
    std::map<int, RefAlt>::iterator posIter = chrVariant[chr].find(posistion);
    // empty position
    if( posIter == chrVariant[chr].end() )
        return false;
        
    return true;
}

void VariantParser::filterSNP(std::string chr, std::vector<ReadVariant> &readVariantVec, std::string &chr_reference){
    
    // pos, <allele, <strand, True>
    std::map<int, std::map<int, std::map<int, bool> > > posAlleleStrand;
    std::map< int, bool > methylation;
    
    // iter all variant, record the strand contained in each SNP
    for( auto readSNPVecIter : readVariantVec ){
        // tag allele on forward or reverse strand
        for(auto variantIter : readSNPVecIter.variantVec ){
            posAlleleStrand[variantIter.position][variantIter.allele][readSNPVecIter.is_reverse] = true;
        }
    }
    // iter all SNP, both alleles that require SNP need to appear in the two strand
    for(auto pos: posAlleleStrand){
        // this position contain two allele, REF allele appear in the two strand, ALT allele appear in the two strand
        if(pos.second.size() == 2 && pos.second[0].size() == 2 && pos.second[1].size() == 2 ){
            // high confident SNP
        }
        else{
            methylation[pos.first] = true;
        }
    }

    // Filter SNPs that are not easy to phasing due to homopolymer
    // get variant list
    std::map<std::string, std::map<int, RefAlt> >::iterator chrIter =  chrVariant.find(chr);
    std::map< int, bool > errorProneSNP;
    
    if( chrIter != chrVariant.end() ){
        std::map<int, int> consecutiveAllele;
        // iter all SNP and tag homopolymer
        for(std::map<int, RefAlt>::iterator posIter = (*chrIter).second.begin(); posIter != (*chrIter).second.end(); posIter++ ){
            consecutiveAllele[(*posIter).first] = homopolymerLength((*posIter).first, chr_reference);
        }
        std::map<int, RefAlt>::iterator currSNPIter = (*chrIter).second.begin();
        std::map<int, RefAlt>::iterator nextSNPIter = std::next(currSNPIter,1);
        // check whether each SNP pair falls in an area that is not easy to phasing
        while( currSNPIter != (*chrIter).second.end() && nextSNPIter != (*chrIter).second.end() ){
            int currPos = (*currSNPIter).first;
            int nextPos = (*nextSNPIter).first;
            // filter one of SNP if this SNP pair falls in homopolymer and distance<=2
            if( consecutiveAllele[currPos] >= 3 && consecutiveAllele[nextPos] >= 3 && std::abs(currPos-nextPos)<=2 ){
                errorProneSNP[nextPos]=true;
                nextSNPIter = (*chrIter).second.erase(nextSNPIter);
                continue;
            }
            
            currSNPIter++;
            nextSNPIter++;
        }
        
    }
    
    
    // iter all reads
    for( std::vector<ReadVariant>::iterator readSNPVecIter = readVariantVec.begin() ; readSNPVecIter != readVariantVec.end() ; readSNPVecIter++ ){
        // iter all SNPs in this read
        for(std::vector<Variant>::iterator variantIter = (*readSNPVecIter).variantVec.begin() ; variantIter != (*readSNPVecIter).variantVec.end() ; ){
            std::map< int, bool >::iterator delSNPIter = methylation.find((*variantIter).position);
            std::map< int, bool >::iterator homoIter = errorProneSNP.find((*variantIter).position);
            
            if( delSNPIter != methylation.end() ){
                variantIter = (*readSNPVecIter).variantVec.erase(variantIter);
            }
            else if( homoIter != errorProneSNP.end() ){
                variantIter = (*readSNPVecIter).variantVec.erase(variantIter);
            }
            else{
                variantIter++;
            }
        }
    }
}

std::map<int, bool> VariantParser::getHomopolymeVariants(std::string chrName){
    std::map<int, bool> targetVariants;
    std::map<std::string, std::map<int, bool> >::iterator chrIter = chrVariantHomopolymer.find(chrName);
    
    if( chrIter != chrVariantHomopolymer.end() )
        targetVariants =  (*chrIter).second;
        
    return targetVariants;
}

int VariantParser::homopolymerEnd(int snp_pos, const std::string &ref_string){
    char element = ref_string.at(snp_pos);
    int ref_len = ref_string.length();
    int pos = snp_pos;
    
    while( ref_string.at(pos) == element ){
        pos++;
        if( pos + 1 >= ref_len )
            return pos;
    }
    
    return pos-1;
}

SVParser::SVParser(std::string variantFile, VariantParser &in_snpFile):variantFile(variantFile){
    
    snpFile = &in_snpFile;
    
    if( variantFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressParser(variantFile);
    }
    else if( variantFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressParser(variantFile);
    }

    // erase SV pos if this pos appear two or more times
    for(std::map<std::string, std::map<int, bool> >::iterator chrIter = posDuplicate.begin(); chrIter != posDuplicate.end() ; chrIter++){
        for(std::map<int, bool>::iterator posIter = (*chrIter).second.begin() ; posIter != (*chrIter).second.end() ; posIter++ ){
            if( (*posIter).second == true){
                std::map<int, std::map<std::string ,bool> >::iterator erasePosIter = chrVariant[(*chrIter).first].find((*posIter).first);
                if(erasePosIter != chrVariant[(*chrIter).first].end()){
                    chrVariant[(*chrIter).first].erase(erasePosIter);
                }
            }
        }
    }
}
SVParser::~SVParser(){ 
}

void SVParser::compressParser(std::string &variantFile){
    gzFile file = gzopen(variantFile.c_str(), "rb");
    if(variantFile=="")
        return;
    if(!file){
        std::cout<< "Fail to open vcf: " << variantFile << "\n";
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
void SVParser::unCompressParser(std::string &variantFile){
    std::ifstream originVcf(variantFile);
    if(variantFile=="")
        return;
    if(!originVcf.is_open()){
        std::cout<< "Fail to open vcf: " << variantFile << "\n";
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
void SVParser::parserProcess(std::string &input){
    if( input.find("#")!= std::string::npos){

    }
    else if( input.find("#")== std::string::npos ){
        std::istringstream iss(input);
        std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

        if( fields.size() == 0 )
            return;
        // trans to 0-base
        int pos = std::stoi( fields[1] ) - 1;
        std::string chr = fields[0];
                
        // find GT flag
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
        bool filter = false;
                
        // homo GT
        if(fields[9][modifu_start]==fields[9][modifu_start+2]){
            filter = true;
        }
        // conflict pos with SNP
        if( (*snpFile).findSNP(chr,pos) ){
            filter = true;
        }
                
        std::map<int, bool>::iterator posIter = posDuplicate[chr].find(pos);
        // conflict pos with SV
        if( posIter == posDuplicate[chr].end() )
            posDuplicate[chr][pos] = false;
        else{
            posDuplicate[chr][pos] = true;
            filter = true;
        }
                
        if(filter){
            return;
        }

        // get read INFO
        int read_pos = fields[7].find("RNAMES=");
        read_pos = fields[7].find("=",read_pos);
        read_pos++;
                
        int next_field = fields[7].find(";",read_pos);
        std::string totalRead = fields[7].substr(read_pos,next_field-read_pos);
        std::stringstream totalReadStream(totalRead);
                
        std::string read;
        while(std::getline(totalReadStream, read, ','))
        {
           chrVariant[chr][pos][read]= true;
        }
    }
}


std::map<int, std::map<std::string ,bool> > SVParser::getVariants(std::string chrName){
    std::map<int, std::map<std::string ,bool> > targetVariants;
    std::map<std::string, std::map<int, std::map<std::string ,bool> > >::iterator chrIter = chrVariant.find(chrName);
    
    if( chrIter != chrVariant.end() )
        targetVariants =  (*chrIter).second;
        
    return targetVariants;
}

void SVParser::writeResult(std::string resultPrefix, PhasingResult phasingResult){
    
    if( variantFile.find("gz") != std::string::npos ){
        // .vcf.gz 
        compressInput(resultPrefix,phasingResult);
    }
    else if( variantFile.find("vcf") != std::string::npos ){
        // .vcf
        unCompressInput(resultPrefix,phasingResult);
    }
    return;
}

void SVParser::compressInput(std::string resultPrefix, PhasingResult phasingResult){
    
    gzFile file = gzopen(variantFile.c_str(), "rb");
    std::ofstream resultVcf(resultPrefix+"_SV.vcf");
    
    if(!resultVcf.is_open()){
        std::cout<< "Fail to open write file: " << resultPrefix+"_SV.vcf" << "\n";
    }
    else if(!file){
        std::cout<< "Fail to open vcf: " << variantFile << "\n";
    }
    else{
        bool ps_def = false;
              
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
                writeLine(input, ps_def, resultVcf, phasingResult);
            }
            // any trailing data in [eol, end) now is a partial line
            offset = std::copy(cur, end, buffer);
        }
        gzclose (file);
    }
}
void SVParser::unCompressInput(std::string resultPrefix, PhasingResult phasingResult){
    std::ifstream originVcf(variantFile);
    std::ofstream resultVcf(resultPrefix+"_SV.vcf");
    
    if(!resultVcf.is_open()){
        std::cout<< "Fail to open write file: " << resultPrefix+"_SV.vcf" << "\n";
    }
    else if(!originVcf.is_open()){
        std::cout<< "Fail to open vcf: " << variantFile << "\n";
    }
    else{
        bool ps_def = false;
        std::string input;
        while(! originVcf.eof() ){
            std::getline(originVcf, input);
            if(input.length() == 0)
                continue;

            writeLine(input, ps_def, resultVcf, phasingResult);
        }
    }
}
void SVParser::writeLine(std::string &input, bool &ps_def, std::ofstream &resultVcf, PhasingResult &phasingResult){
    // header
    if( input.find("#")!= std::string::npos){
        // avoid double definition
        if( input.find("ID=PS")!= std::string::npos){
            ps_def = true;
        }
        // format line 
        if(input.find("##")== std::string::npos && ps_def == false){
            resultVcf <<  "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">\n";
        }
        resultVcf << input << "\n";
    }
    else if( input.find("#")== std::string::npos ){
        std::istringstream iss(input);
        std::vector<std::string> fields((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());

        std::string chr = fields[0];

        if( fields.size() == 0 )
            return;

        int pos = std::stoi( fields[1] );
        int posIdx = pos - 1 ;
        std::string key = fields[0] + "_" + std::to_string(posIdx);

        PhasingResult::iterator psElementIter =  phasingResult.find(key);
        std::map<int, std::map<std::string ,bool> >::iterator posIter = chrVariant[chr].find(posIdx);
                
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
            int modifu_start = 0;
            for(unsigned int i =0; i < fields[9].length() ; i++){
                if( current_colon >= colon_pos )
                    break;
                if(fields[9][i]==':')
                    current_colon++;  
                modifu_start++;
            }
            if(fields[9][modifu_start+1] == '|'){
                // direct modify GT value
                if(fields[9][modifu_start] > fields[9][modifu_start+2]){
                    fields[9][modifu_start+1] = fields[9][modifu_start];
                    fields[9][modifu_start] = fields[9][modifu_start+2];
                    fields[9][modifu_start+2] =fields[9][modifu_start+1];
                }
                fields[9][modifu_start+1] = '/';
            }        
        }
                
        // this pos is phase and exist in map
        if( psElementIter != phasingResult.end() && posIter != chrVariant[chr].end() ){
                    
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
            int modifu_start = 0;
            for(unsigned int i =0; i < fields[9].length() ; i++){
                if( current_colon >= colon_pos )
                    break;
                if(fields[9][i]==':')
                    current_colon++;  
                modifu_start++;
            }
            // direct modify GT value
            fields[9][modifu_start] = (*psElementIter).second.RAstatus[0];
            fields[9][modifu_start+1] = '|';
            fields[9][modifu_start+2] = (*psElementIter).second.RAstatus[2];
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

BamParser::BamParser(std::string inputChrName, std::string inputBamFile, VariantParser snpMap, SVParser svFile):chrName(inputChrName),BamFile(inputBamFile){
    // use chromosome to find recorded snp map
    currentVariants = snpMap.getVariants(chrName);
    // set skip variant start iterator
    firstVariantIter = currentVariants.begin();
    if( firstVariantIter == currentVariants.end() ){
        std::cerr<< "error chromosome name or empty map.\n";
        exit(1);
    }
    else{
        
    }
    currentSV = svFile.getVariants(chrName);
    firstSVIter = currentSV.begin();
}
BamParser::~BamParser(){
}

void BamParser::direct_detect_alleles(int lastSNPPos, PhasingParameters params, std::vector<ReadVariant> &readVariantVec, const std::string &ref_string){
    int numThreads = params.numThreads;
    char *bamList[numThreads];
    int chunkSize = lastSNPPos/numThreads + numThreads;
    
    // create list. i.e. chr1:100-200
    for(int i=0;i<numThreads;i++){
        std::string tmp = chrName + ":" + std::to_string(i*chunkSize+1) + "-" + std::to_string((i+1)*chunkSize);
        bamList[i] = new char[tmp.length()+1];
        strcpy(bamList[i], tmp.c_str());
    }
    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};
    // open bam file
    samFile *fp_in = hts_open(BamFile.c_str(),"r"); 
    // read header
	bam_hdr_t *bamHdr = sam_hdr_read(fp_in); 
    // initialize an alignment
	bam1_t *aln = bam_init1(); 
    hts_idx_t *idx = NULL;
    
    if ((idx = sam_index_load(fp_in, BamFile.c_str())) == 0) {
        std::cout<<"ERROR: Cannot open index for bam file\n";
        exit(1);
    }
    
    hts_itr_multi_t *iter;

    if( (iter = sam_itr_regarray(idx, bamHdr, bamList, numThreads)) == 0){
        std::cout<<"ERROR: Cannot open iterator for region %s for bam file\n";
        exit(1);
    }
    
    int result;
    
    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
    }
    hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &threadPool);
    
    while ((result = sam_itr_multi_next(fp_in, iter, aln)) >= 0) { 
        int flag = aln->core.flag;

        if (    aln->core.qual < 20  // mapping quality
             || (flag & 0x4)   != 0  // read unmapped
             || (flag & 0x100) != 0  // secondary alignment. repeat. 
                                     // A secondary alignment occurs when a given read could align reasonably well to more than one place.
             || (flag & 0x400) != 0  // duplicate 
             // (flag & 0x800) != 0  // supplementary alignment
                                     // A chimeric alignment is represented as a set of linear alignments that do not have large overlaps.
             ){ 
            std::string str = bamHdr->target_name[aln->core.tid];
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
        
        get_snp(tmp, readVariantVec, ref_string, params.isONT);
        free(tmp.quality);
        free(tmp.qseq);
        free(tmp.op);
        free(tmp.ol);
	}
	hts_idx_destroy(idx);
    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
	sam_close(fp_in);
    hts_tpool_destroy(threadPool.pool);
}

bool BamParser::continueHomopolimer(int snp_pos, const std::string &ref_string){
    int homopolymer_length = 1;
    char element = ref_string.at(snp_pos);
    bool neighbor_is_homopolymer = false;
    
    int pos = snp_pos-1;
    while( ref_string.at(pos) == element ){
        pos--;
        homopolymer_length++;
    }
    
    if( homopolymerLength(pos,ref_string) >= 3 )
        neighbor_is_homopolymer = true;
        
    pos = snp_pos+1;
    while( ref_string.at(pos) == element ){
        pos++;
        homopolymer_length++;
    }
    
    if( homopolymerLength(pos,ref_string) >= 3 )
        neighbor_is_homopolymer = true;
    
    if( homopolymer_length >=3 && neighbor_is_homopolymer)
        return true;
    
    return false;
}

void BamParser::get_snp(Alignment align, std::vector<ReadVariant> &readVariantVec,const std::string &ref_string, bool isONT){

    ReadVariant *tmpReadResult = new ReadVariant();
        (*tmpReadResult).read_name = align.qname;
        (*tmpReadResult).source_id = align.chr;
        (*tmpReadResult).reference_start = align.refStart;
        (*tmpReadResult).is_reverse = align.is_reverse;
    
        // Skip variants that are to the left of this read
        while( firstVariantIter != currentVariants.end() && (*firstVariantIter).first < align.refStart )
            firstVariantIter++;
        
        // Skip structure variants that are to the left of this read
        while( firstSVIter != currentSV.end() && (*firstSVIter).first < align.refStart )
            firstSVIter++;
        
        // position relative to reference
        int ref_pos = align.refStart;
        // position relative to read
        int query_pos = 0;
        // translation char* to string;
        std::string qseq = align.qseq;
        // translation char* to string;
        
        // set variant start for current alignment
        std::map<int, RefAlt>::iterator currentVariantIter = firstVariantIter;
        
        // set structure variant start for current alignment
         std::map<int, std::map<std::string ,bool> >::iterator currentSVIter = firstSVIter;
        
        std::map<int, RefAlt>::reverse_iterator last = currentVariants.rbegin();
        
        // reading cigar to detect snp on this read
        for(int i = 0; i < align.cigar_len ; i++ ){
            int cigar_op = align.op[i];
            int length   = align.ol[i];
            
            // iterator next variant
            while( currentVariantIter != currentVariants.end() && (*currentVariantIter).first < ref_pos ){
                currentVariantIter++;
            }
           
            while( currentSVIter != currentSV.end() && (*currentSVIter).first < ref_pos )
                currentSVIter++;
            
            if( currentSVIter != currentSV.end() ){
                if( (*currentSVIter).first < (*currentVariantIter).first ){
                    std::map<std::string ,bool>::iterator readIter = currentSV[(*currentSVIter).first].find(align.qname);
                    // If this read not contain SV, it means this read is the same as reference genome.
                    // default this read the same as ref
                    int allele = 0; 
                    // this read contain SV.
                    if( readIter != currentSV[(*currentSVIter).first].end() )
                        allele = 1;
                    
                    // generate SV read count
                    if(allele==0){
                        std::map<int, int>::iterator svCountIter = noSV.find((*currentSVIter).first);
                        if( svCountIter == noSV.end() ){
                            noSV[(*currentSVIter).first] = 0;
                        }
                        else{
                            noSV[(*currentSVIter).first]++;
                        }
                    }
                    else{
                        std::map<int, int>::iterator svCountIter = occurSV.find((*currentSVIter).first);
                        if( svCountIter == occurSV.end() ){
                            occurSV[(*currentSVIter).first] = 0;
                        }
                        else{
                            occurSV[(*currentSVIter).first]++;
                        }
                    }
                    
                    // push this SV to vector
                    Variant *tmpVariant = new Variant((*currentSVIter).first, allele, -1 );
                    (*tmpReadResult).variantVec.push_back( (*tmpVariant) );
                    // next SV iter     
                    currentSVIter++;
                }
            }
            // CIGAR operators: MIDNSHP=X correspond 012345678
            // 0: alignment match (can be a sequence match or mismatch)
            // 7: sequence match
            // 8: sequence mismatch
            if( cigar_op == 0 || cigar_op == 7 || cigar_op == 8 ){
                
                while( currentVariantIter != currentVariants.end() && (*currentVariantIter).first < ref_pos + length){

                    int refAlleleLen = (*currentVariantIter).second.Ref.length();
                    int altAlleleLen = (*currentVariantIter).second.Alt.length();
                    int offset = (*currentVariantIter).first - ref_pos;

                    std::string base = qseq.substr(query_pos + offset, 1);
                    
                    int allele = -1;
  
                    if( refAlleleLen == 1 && altAlleleLen == 1){
                        if( base == (*currentVariantIter).second.Ref )
                            allele = 0;
                        else if( base == (*currentVariantIter).second.Alt )
                            allele = 1;
                    }
                    else{
                        // skip this position
                        //std::cout<< "not snp: " <<(*currentVariantIter).first + 1 << "\n";
                        //exit(1);
                    }            
                    
                    if( allele != -1 ){
                        int base_q = align.quality[query_pos + offset];
                        // record snp result
                        Variant *tmpVariant = new Variant((*currentVariantIter).first, allele, base_q);
                        (*tmpReadResult).variantVec.push_back( (*tmpVariant) );
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
                
                // If a reference is given
                // it will determine whether the SNP falls in the homopolymer
                // and start the processing of the SNP fall in the alignment GAP
                
                if(ref_string != "" && isONT){
                    int del_len = length;
                    if ( ref_pos + del_len + 1 == (*currentVariantIter).first ){
                        if( homopolymerLength((*currentVariantIter).first , ref_string) >=3 ){
                            // special case
                        }
                    }
                    else if( (*currentVariantIter).first >= ref_pos  && (*currentVariantIter).first < ref_pos + del_len ){
                        // check snp in homopolymer
                        if( homopolymerLength((*currentVariantIter).first , ref_string) >=3 ){
                            
                            int refAlleleLen = (*currentVariantIter).second.Ref.length();
                            int altAlleleLen = (*currentVariantIter).second.Alt.length();
                            
                            // get the next match
                            std::string base = qseq.substr(query_pos , 1);

                            int allele = -1;

                            if( refAlleleLen == 1 && altAlleleLen == 1){
                                if( base == (*currentVariantIter).second.Ref )
                                    allele = 0;
                                else if( base == (*currentVariantIter).second.Alt )
                                    allele = 1;
                            }
                            if(allele != -1){
                                int base_q = align.quality[query_pos];
                                Variant *tmpVariant = new Variant((*currentVariantIter).first, allele, base_q/2);
                                (*tmpReadResult).variantVec.push_back( (*tmpVariant) );
                                currentVariantIter++;
                            }
                            
                        }
                        else{
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
                std::cerr<< "alignment find unsupported CIGAR operation from read: " << align.qname << "\n";
                exit(1);
            }

        }
        
        if( (*tmpReadResult).variantVec.size() > 0 )
            readVariantVec.push_back((*tmpReadResult));
}

void BamParser::svFilter(std::vector<ReadVariant> &readVariantVec){
    
    std::map< int, bool > lowConfidentSV;
    
    for(std::map<int, std::map<std::string ,bool> >::iterator svIter = currentSV.begin() ; svIter != currentSV.end() ; svIter++ ){
        int currSVPos= (*svIter).first;
        
        std::map<int, int>::iterator occurSVIter = occurSV.find(currSVPos);
        std::map<int, int>::iterator noSVIter = noSV.find(currSVPos);
        
        double occurSVNum = occurSVIter != occurSV.end() ? occurSV[currSVPos] : 0;
        double noSVNum    = noSVIter    != noSV.end()    ? noSV[currSVPos] : 0;
        
        //std::cout<<"SV\t"<< currSVPos << "\t"<< occurSVNum << "\t" << noSVNum << "\t" << occurSVNum+noSVNum << "\t" << std::min(occurSVNum,noSVNum)/(occurSVNum+noSVNum) << "\n";
        
        if( std::min(occurSVNum,noSVNum)/(occurSVNum+noSVNum)<=0.1 ){
        //if( occurSVNum + noSVNum < 10 || occurSVNum == 0 || noSVNum == 0 ){
            lowConfidentSV[currSVPos]=true;
        }
    }
    
    if( currentSV.size() > 0)
    {
        // iter all reads
        for( std::vector<ReadVariant>::iterator readSNPVecIter = readVariantVec.begin() ; readSNPVecIter != readVariantVec.end() ; readSNPVecIter++ ){
            // iter all SNPs in this read
            for(std::vector<Variant>::iterator variantIter = (*readSNPVecIter).variantVec.begin() ; variantIter != (*readSNPVecIter).variantVec.end() ; ){
                std::map< int, bool >::iterator svIter = lowConfidentSV.find((*variantIter).position);
                
                if( svIter != lowConfidentSV.end() ){
                    variantIter = (*readSNPVecIter).variantVec.erase(variantIter);
                    //std::cout<<"del\t" << (*svIter).first << "\n";
                }
                else{
                    variantIter++;
                }
            }
        }
    }
}





