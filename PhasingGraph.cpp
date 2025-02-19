#include "PhasingGraph.h"
//SubEdge

SubEdge::SubEdge():readCount(0){ 
    refRead = new std::map<int, std::vector<std::string> >;
    altRead = new std::map<int, std::vector<std::string> >;
    refQuality = new std::map<int, int>;
    altQuality = new std::map<int, int>;
    refReadCount = new std::map<int, float>;
    altReadCount = new std::map<int, float>;
}

SubEdge::~SubEdge(){ 
}

void SubEdge::destroy(){
    delete refRead;
    delete altRead;
    delete refQuality;
    delete altQuality;
    delete refReadCount;
    delete altReadCount;
}

void SubEdge::addSubEdge(int currentQuality, Variant connectNode, std::string readName, int baseQuality, double edgeWeight){
    // target noded is REF allele
    if(connectNode.allele == 0 ){
        // debug, this parameter will record the names of all reads between two points
        //(*refRead)[connectNode.position].push_back(readName);
        /*// quality sum
        std::map<int, int>::iterator rqIter = refQuality->find(connectNode.position);
        if( rqIter == refQuality->end() ){
            (*refQuality)[connectNode.position] = currentQuality + connectNode.quality;
        }
        else{
            (*refQuality)[connectNode.position] += currentQuality + connectNode.quality;
        }*/

	//if the base quality on both snps is high enough, the edge has normal weight
        if ( currentQuality >= baseQuality && connectNode.quality >= baseQuality)
            (*refReadCount)[connectNode.position]++;
        else
            (*refReadCount)[connectNode.position] = (*refReadCount)[connectNode.position] + edgeWeight ;
        
        
	//(*refReadCount)[connectNode.position]++;
    }
    // target noded is ALT allele
    else if(connectNode.allele == 1 ){
        // debug, this parameter will record the names of all reads between two points
        // (*altRead)[connectNode.position].push_back(readName);
        /*// quality sum
        std::map<int, int>::iterator aqIter = altQuality->find(connectNode.position);
        if( aqIter == altQuality->end() ){
            (*altQuality)[connectNode.position] = currentQuality + connectNode.quality;
        }
        else{
            (*altQuality)[connectNode.position] += currentQuality + connectNode.quality;
        }*/
       
	//if the base quality on both snps is high enough, the edge has normal weight 
        if ( currentQuality >= baseQuality && connectNode.quality >= baseQuality)
            (*altReadCount)[connectNode.position]++;
        else
            (*altReadCount)[connectNode.position] = (*altReadCount)[connectNode.position] + edgeWeight ;
        
	//(*altReadCount)[connectNode.position]++;
    }
    readCount++;
}

std::pair<float,float> SubEdge::BestPair(int targetPos){
    return std::make_pair( getRefReadCount(targetPos), getAltReadCount(targetPos) );
}

float SubEdge::getRefReadCount(int targetPos){
    std::map<int, float>::iterator posIter = refReadCount->find(targetPos);
    if( posIter != refReadCount->end() ){
        return (*refReadCount)[targetPos];
    }
    return 0;
}

float SubEdge::getAltReadCount(int targetPos){
    std::map<int, float>::iterator posIter = altReadCount->find(targetPos);
    if( posIter != altReadCount->end() ){
        return (*altReadCount)[targetPos];
    }
    return 0;
}

std::vector<std::string> SubEdge::showEdge(std::string message){
    std::vector<std::string> result;
    for(std::map<int, float >::iterator edgeIter = refReadCount->begin() ; edgeIter != refReadCount->end() ; edgeIter++ ){
        result.push_back(message +" -> ref_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second) + "];");
    }
    for(std::map<int, float >::iterator edgeIter = altReadCount->begin() ; edgeIter != altReadCount->end() ; edgeIter++ ){
        result.push_back(message +" -> alt_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second) + "];");
    }
    return result;
}

std::vector<std::pair<int,int>> SubEdge::getConnectPos(){
    std::vector<std::pair<int,int>> result;
    for(std::map<int, float >::iterator edgeIter = refReadCount->begin() ; edgeIter != refReadCount->end() ; edgeIter++ ){
        result.push_back( std::make_pair( (*edgeIter).first, 0 ) );
    }
    for(std::map<int, float >::iterator edgeIter = altReadCount->begin() ; edgeIter != altReadCount->end() ; edgeIter++ ){
        result.push_back( std::make_pair( (*edgeIter).first, 1 ) );
    }
    return result;
}

int SubEdge::getQuality(PosAllele targetPos){
    // target is Ref allele
    if( targetPos.second == 1 ){
        std::map<int, int>::iterator qIter = refQuality->find(targetPos.first);
        if( qIter == refQuality->end() )
            return 0;
        else
            return (*refQuality)[targetPos.first];
    }
    // target is Alt allele
    if( targetPos.second == 2 ){
        std::map<int, int>::iterator qIter = altQuality->find(targetPos.first);
        if( qIter == altQuality->end() )
            return 0;
        else
            return (*altQuality)[targetPos.first];
    }
    return 0;
}

int SubEdge::getAvgQuality(PosAllele targetPos){
    // target is Ref allele
    if( targetPos.second == 1 ){
        std::map<int, int>::iterator qIter = refQuality->find(targetPos.first);
        if( qIter == refQuality->end() )
            return 0;
        else
            return (*refQuality)[targetPos.first]/(*refReadCount)[targetPos.first];
    }
    // target is Alt allele
    if( targetPos.second == 2 ){
        std::map<int, int>::iterator qIter = altQuality->find(targetPos.first);
        if( qIter == altQuality->end() )
            return 0;
        else
            return (*altQuality)[targetPos.first]/(*altReadCount)[targetPos.first];
    }
    return 0;
}

VoteResult::VoteResult( int currPos, float variantweight ) {
    Pos = currPos ;
    weight = variantweight ;
}

VariantEdge::VariantEdge(int inCurrPos){
    currPos = inCurrPos;
    alt = new SubEdge();
    ref = new SubEdge();
}

//VariantEdge
std::pair<PosAllele,PosAllele> VariantEdge::findBestEdgePair(int targetPos, bool isONT, double edgeThreshold, bool debug, std::map<int,int> &variantType, VoteResult &vote){
    std::pair<float,float> refBestPair  = ref->BestPair(targetPos);
    std::pair<float,float> altBestPair  = alt->BestPair(targetPos);
    // get the weight of each pair
    float rr = refBestPair.first;
    float ra = refBestPair.second;
    float ar = altBestPair.first;
    float aa = altBestPair.second;
    // initialize the edge connection
    // -1 : not connect
    int refAllele = -1;
    int altAllele = -1;
    
    double edgeSimilarRatio = (double)std::min((rr+aa),(ar+ra)) / (double)std::max((rr+aa),(ar+ra));
 
    if( rr + aa > ra + ar ){
        // RR conect
        refAllele = 1;
        altAllele = 2;
    }
    else if( rr + aa < ra + ar ){
        // RA connect
        refAllele = 2;
        altAllele = 1;
    }
    else if( rr + aa == ra + ar ){
        // no connect 
        // not sure which is better
    }

    //VarintType < 0=SNP 1=SV 2=MOD 3=INDEL 4=tandem repeat INDEL >
    if((variantType[currPos] == 0 && variantType[targetPos] == 2)||(variantType[currPos] == 2 && variantType[targetPos] == 0)){
        edgeThreshold = 0.3;
        if((rr+ra+ar+aa) < 1){
            edgeThreshold = -1;
        }
    }

    if( edgeSimilarRatio > edgeThreshold ){
        refAllele = -1;
        altAllele = -1;
    }


    if(debug){
        std::cout << currPos + 1 << "\t->\t" << targetPos + 1 << "\t|rr aa | ra ar\t" << "\t" << rr << "\t" << aa << "\t" << ra << "\t" << ar  << "\n";
    }

    // the lower the edgeSimilarRatio means the higher reads consistency, and we will make the weight bigger if the reads consistency is high enough
    else if ( (edgeSimilarRatio <= 0.1 && (rr + aa + ra + ar) >= 1)  || ((rr+aa)<1&&(ra+ar)>=1) || ((rr+aa)>=1&&(ra+ar)<1) ) {
        vote.weight = 20 ;
    }

    vote.para = rr + aa ;
    vote.cross = ra + ar ;
    vote.ESR = edgeSimilarRatio ;

    // create edge pairs
    PosAllele refEdge = std::make_pair( targetPos, refAllele );
    PosAllele altEdge = std::make_pair( targetPos, altAllele );
    // return edge pair
    return std::make_pair( refEdge, altEdge );
}

std::pair<float,float> VariantEdge::findNumberOfRead(int targetPos){
    std::pair<float,float> refBestPair  = ref->BestPair(targetPos);
    std::pair<float,float> altBestPair  = alt->BestPair(targetPos);
    // get the weight of each pair
    float rr = refBestPair.first;
    float ra = refBestPair.second;
    float ar = altBestPair.first;
    float aa = altBestPair.second;
    return std::make_pair( rr + aa , ra +ar );
}

//BlockRead
void BlockRead::recordRead(std::string readName){
    std::map<std::string,int>::iterator readIter = readVec.find(readName);
    if( readIter == readVec.end() )
        readVec[readName] = 1;
    else
        readVec[readName]++;
}

//Handle the special case which One Long Read provides wrong info repeatedly
std::pair<float,float> VairiantGraph::Onelongcase( std::vector<VoteResult> vote ){

    int counter = 0 ;
    float h1 = 0 ;
    float h2 = 0 ;

    // iterate all the voting that previous variants provide
    for (std::vector<VoteResult>::size_type i = 0 ; i < vote.size() ; i++ ) {

	// count the votes that refer to only one read
        if ( (vote[i].para+vote[i].cross) <= 1 ) {
            counter++ ;
        }
	// we will only count the votes that is not INDEL and have lower ESR beacause the INDEL is the variant has higher error rate and the lower ESR means higher reads consistency,
        else if ( vote[i].ESR < 0.2 && vote[i].weight >= 1 && (*variantType)[vote[i].Pos] != 3 ) {
            if ( vote[i].hap == 1 ) {
                h1+=vote[i].weight ;
            }
            else if ( vote[i].hap == 2 ) {
                h2+=vote[i].weight ;
            }
        }
    }

    //if there has less than three variants use one read to vote we cancel the mechanism
    if ( counter <= 3 || (h1==0&&h2==0) ) {
        return std::make_pair( -1 , -1 ) ;
    }
    else {
        return std::make_pair( h1 , h2 ) ;
    }

}

//VairiantGraph
void VairiantGraph::edgeConnectResult(){
    //current variant position, haplotype (1 or 2), previous variants' voting information 
    std::map<int, std::vector<VoteResult> > *hpCountMap3 = new std::map<int, std::vector<VoteResult> > ;
    // current variant position, haplotype (1 or 2), previous variants' voting result
    std::map<int, std::map<int,float> > *hpCountMap2 = new std::map<int, std::map<int,float> > ;
    // current snp, haplotype (1 or 2), support snp
    std::map<int, std::map<int,std::vector<int> > > *hpCountMap = new std::map<int, std::map<int,std::vector<int> > >;
    // current snp, result haplotype (1 or 2)
    std::map<int,int> *hpResult = new std::map<int,int>;
    // < block start, <snp in this block> >
    std::map<int,std::vector<int> > *phasedBlocks = new std::map<int,std::vector<int> >;
    
    int blockStart = -1;
    int currPos = -1;
    int nextPos = -1;
    int lastConnectPos = -1;

    // Visit all position and assign SNPs to haplotype.
    // Avoid recording duplicate information,
    // only one of the two alleles needs to be used for each SNP
    for(std::map<int,ReadBaseMap*>::iterator variantIter = totalVariantInfo->begin() ; variantIter != totalVariantInfo->end() ; variantIter++ ){
        // check next position
        std::map<int,ReadBaseMap*>::iterator nextNodeIter = std::next(variantIter, 1);
        if( nextNodeIter == totalVariantInfo->end() ){
             break;
        }
        
        currPos = variantIter->first;
        nextPos = nextNodeIter->first;
        
        // There should not be a large distance between any two variants, 
        // with the default being a distance of 300000bp, equivalent to one centromere length.
        if(std::abs(nextPos-currPos) > params->distance ){
            continue;
        }
        
        // get the number of HP1 and HP2 supported reference allele
        //int h1 = (*hpCountMap)[currPos][1].size();
        //int h2 = (*hpCountMap)[currPos][2].size();
        float h1 = (*hpCountMap2)[currPos][1] ;
	    float h2 = (*hpCountMap2)[currPos][2] ;

        //std::cout << currPos+1 << "\t" << "h1:" << "\t" << h1 << "\t" << "h2:" << "\t" << h2 << "\n" ;

	//Handle the special case which One Long Read provides wrong info repeatedly
        std::pair<float, float> special = Onelongcase( (*hpCountMap3)[currPos] ) ;
	    if ( special.first != -1 ) {
           h1 = special.first ;
           h2 = special.second ;
        }

        // new block, set this position as block start 
        if( h1 == h2 ){
            // No new blocks should be created if the next SNP has already been picked up
            if( currPos < lastConnectPos ){
                continue;
            }
            
            blockStart = currPos;
            (*phasedBlocks)[blockStart].push_back(currPos);
            (*hpResult)[currPos] = 1;
        }
        else{
            int currHP = ( h1 > h2 ? 1 : 2 );
            (*hpResult)[currPos] = currHP;
            (*phasedBlocks)[blockStart].push_back(currPos);
        }
        // Check if there is no edge from current node
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList->find( currPos );
        if( edgeIter==edgeList->end() ){
            continue;
        }
        
        // check connect between surrent SNP and next n SNPs
        for(int i = 0 ; i < params->connectAdjacent ; i++ ){
	    VoteResult vote(currPos, 1); //used to store previous 20 variants' voting information

        // consider reads from the currnt SNP and the next (i+1)'s SNP
        std::pair<PosAllele,PosAllele> tmp = edgeIter->second->findBestEdgePair(nextNodeIter->first, params->isONT, params->edgeThreshold, false, *variantType, vote);

	    // if the target is a danger indel change its weight to 0.1
	    if ( (*variantType)[currPos] == 4 ) {
                vote.weight = 0.1 ;
            }

            // -1 : no connect  
            //  1 : the haplotype of next (i+1)'s SNP are same as previous
            //  2 : the haplotype of next (i+1)'s SNP are different as previous
            if( tmp.first.second != -1 ){
                // record the haplotype resut of next (i+1)'s SNP
                if( (*hpResult)[currPos] == 1 ){
                    if( tmp.first.second == 1 ){
                        (*hpCountMap)[nextNodeIter->first][1].push_back(currPos);
		                (*hpCountMap2)[nextNodeIter->first][1] += vote.weight;
                        vote.hap = 1 ;
                    }
                    if( tmp.first.second == 2 ){
                        (*hpCountMap)[nextNodeIter->first][2].push_back(currPos);
		                (*hpCountMap2)[nextNodeIter->first][2] += vote.weight;
		                vote.hap = 2 ;
                    }
                }
                if( (*hpResult)[currPos]==2 ){
                    if( tmp.first.second == 1 ){
                        (*hpCountMap)[nextNodeIter->first][2].push_back(currPos);
                        (*hpCountMap2)[nextNodeIter->first][2] += vote.weight;
                        vote.hap = 2 ;
                    }
                    if( tmp.first.second == 2 ){
                        (*hpCountMap)[nextNodeIter->first][1].push_back(currPos);
                        (*hpCountMap2)[nextNodeIter->first][1] += vote.weight;
                        vote.hap = 1 ;
                    }
                }

                (*hpCountMap3)[nextNodeIter->first].push_back( vote );

                if( params->generateDot ){
                    std::string e1 = std::to_string(currPos+1) + ".1\t->\t" + std::to_string(tmp.first.first+1) + "." + std::to_string(tmp.first.second);
                    std::string e2 = std::to_string(currPos+1) + ".2\t->\t" + std::to_string(tmp.second.first+1) + "." + std::to_string(tmp.second.second);

                    dotResult.push_back(e1);
                    dotResult.push_back(e2);
                }
                
                lastConnectPos = nextNodeIter->first;
            }
            nextNodeIter++;
            if( nextNodeIter == totalVariantInfo->end() ){
                break;
            }
        }
    }

    // outFile.close();
    // loop all block and construct graph
    // Record the phase set(PS) for each variant on the graph and record the haplotype to each variant's allele belongs.
    for(auto blockIter = phasedBlocks->begin() ; blockIter != phasedBlocks->end() ; blockIter++ ){
        // check block size to skip one node island
        if( (*blockIter).second.size()<=1 ){
            continue;
        }
        
        // loop block's node
        // store phasing results include PS and HP
        for(auto currIter = (*blockIter).second.begin() ; currIter != (*blockIter).second.end() ; currIter++ ){
            // check next node
            auto nextIter = std::next(currIter,1);
            if( nextIter == (*blockIter).second.end() ){
                continue;
            }

            PosAllele refStart = std::make_pair((*currIter), 1);
            PosAllele altStart = std::make_pair((*currIter), 2);
            PosAllele refEnd = std::make_pair((*nextIter), 1);
            PosAllele altEnd = std::make_pair((*nextIter), 2);
            
            // store PS results
            (*bkResult)[refStart] = (*blockIter).first + 1;
            (*bkResult)[refEnd]   = (*blockIter).first + 1;
            (*bkResult)[altStart] = (*blockIter).first + 1;
            (*bkResult)[altEnd]   = (*blockIter).first + 1;
            
            // store HP results
            if( currIter == (*blockIter).second.begin() ){
                (*subNodeHP)[refStart] = 0;
                (*subNodeHP)[altStart] = 1;
            }
            
            if( (*hpResult)[(*currIter)] == 0 || (*hpResult)[(*nextIter)] == 0 ){
                
            }
            else if( (*hpResult)[(*currIter)] == (*hpResult)[(*nextIter)] ){
                (*subNodeHP)[refEnd] = (*subNodeHP)[refStart];
                (*subNodeHP)[altEnd] = (*subNodeHP)[altStart];
            }
            else{
                (*subNodeHP)[refEnd] = (*subNodeHP)[altStart];
                (*subNodeHP)[altEnd] = (*subNodeHP)[refStart];
            }
        }
    }
    
    delete hpCountMap;
    delete hpCountMap2;
    delete hpCountMap3;
    delete hpResult;
    delete phasedBlocks;
}

VairiantGraph::VairiantGraph(std::string &in_ref, PhasingParameters &in_params){
    params=&in_params;
    ref=&in_ref;
    totalVariantInfo = new std::map<int,ReadBaseMap*>;
    edgeList = new std::map<int,VariantEdge*>;
    bkResult = new std::map<PosAllele,int>;
    subNodeHP = new std::map<PosAllele,int>;
    variantType = new std::map<int,int>;
    readHpMap = new std::map<std::string,int>;
}

VairiantGraph::~VairiantGraph(){
}

void VairiantGraph::destroy(){
    dotResult.clear();
    dotResult.shrink_to_fit();

    for( auto edgeIter = edgeList->begin() ; edgeIter != edgeList->end() ; edgeIter++ ){
        edgeIter->second->ref->destroy();
        edgeIter->second->alt->destroy();
        delete edgeIter->second->ref;
        delete edgeIter->second->alt;
    }
    
    for( auto variantIter = totalVariantInfo->begin() ; variantIter != totalVariantInfo->end() ; variantIter++ ){
        delete variantIter->second;
    }
    
    delete totalVariantInfo;
    delete edgeList;
    delete bkResult;
    delete subNodeHP;
    delete variantType;
    delete readHpMap;
}

//check if the position is in the range of the cnv
bool VairiantGraph::isPositionInRange(int position, int start, int end){
    return position >= start && position <= end;
}

//calculate the mismatch rate of the cnv
void VairiantGraph::calculateCnvMismatchRate(std::vector<ReadVariant>& in_readVariant, Clip &clip){
    //if the read variant is empty or the cnv is empty, return
    if(in_readVariant.empty() || clip.cnvVec.empty()){
        return;
    }
    size_t cnvIndex = 0;

    for(auto& read : in_readVariant){
        //if the read has no variant, continue
        if(read.variantVec.empty()){
            continue;
        }

        int readStart = read.variantVec.front().position;
        int readEnd = read.variantVec.back().position;

        //find the index of the cnv vector  
        while(cnvIndex > 0 && clip.cnvVec[cnvIndex].first > readStart){
            cnvIndex--;
        }
        size_t i = cnvIndex;

        //iterate through the cnv vector
        while(i < clip.cnvVec.size() && clip.cnvVec[i].first <= readEnd){
            for(const auto& variant : read.variantVec){
                //if the variant position is greater than the end of the cnv, break
                if(variant.position > clip.cnvVec[i].second){
                    break;
                }
                //if the variant position is in the range of the cnv and the allele is reference, increment the mismatch rate
                if(isPositionInRange(variant.position, clip.cnvVec[i].first, clip.cnvVec[i].second) && variant.allele == 1){
                    read.cnv_mmrate_map[clip.cnvVec[i].first]++;
                }
            }
            i++;
        }
        //update the index of the cnv vector
        cnvIndex = i > 0 ? i - 1 : 0;
    }
    
}

//aggregate the mismatch rate of the cnv
void VairiantGraph::aggregateCnvReadMismatchRate(const std::vector<ReadVariant>& in_readVariant, const Clip &clip, std::map<int, std::map<int, std::vector<int>>>& cnvReadMmrate) { 
    //if the read variant is empty or the cnv is empty, return
    if(in_readVariant.empty() || clip.cnvVec.empty()){
        return;
    }

    size_t cnvIndex = 0;

    for(const auto& read : in_readVariant){
        //if the read has no variant, continue
        if(read.variantVec.empty()){
            continue;
        }

        int readStart = read.variantVec.front().position;
        int readEnd = read.variantVec.back().position;

        //find the index of the cnv vector
        while(cnvIndex > 0 && clip.cnvVec[cnvIndex].first > readStart){
            cnvIndex--;
        }

        size_t i = cnvIndex;
        
        //iterate through the cnv vector
        while(i < clip.cnvVec.size() && clip.cnvVec[i].first <= readEnd){
            for(const auto& variant : read.variantVec){
                //if the variant position is greater than the end of the cnv, break
                if(variant.position > clip.cnvVec[i].second){
                    break;
                }
                //if the variant position is in the range of the cnv, push the mismatch rate to the vector to the direct position and allele
                if(isPositionInRange(variant.position, clip.cnvVec[i].first, clip.cnvVec[i].second) && read.cnv_mmrate_map.find(clip.cnvVec[i].first) != read.cnv_mmrate_map.end()){
                    cnvReadMmrate[variant.position][variant.allele].push_back(read.cnv_mmrate_map.at(clip.cnvVec[i].first));
                }
            }
            i++;
        }
        cnvIndex = i > 0 ? i - 1 : 0;
    }
}

//calculate the average mismatch rate of the cnv
void VairiantGraph::calculateAverageMismatchRate(const Clip& clip, const std::map<int, std::map<int, std::vector<int>>>& cnvReadMmrate, std::map<int, double>& missRateMap){

    if(cnvReadMmrate.empty() || clip.cnvVec.empty()){
        return;
    }

    size_t cnvIndex = 0;

    for(const auto& variant : cnvReadMmrate){
        while(cnvIndex > 0 && clip.cnvVec[cnvIndex].first > variant.first){
            cnvIndex--;
        }
        
        size_t i = cnvIndex;
        //iterate through the cnv vector
        while(i < clip.cnvVec.size()){
            //if the variant position is greater than the end of the cnv, break
            if(clip.cnvVec[i].first > variant.first){
                break;
            }
            //if the variant position is in the range of the cnv, calculate the average mismatch rate
            if(isPositionInRange(variant.first, clip.cnvVec[i].first, clip.cnvVec[i].second)){
                auto ref_iter = variant.second.find(0);
                auto alt_iter = variant.second.find(1);
                if(ref_iter != variant.second.end() && alt_iter != variant.second.end()){
                    double AvgRefCnvReadMiss = calculateMean(ref_iter->second);
                    double AvgAltCnvReadMiss = calculateMean(alt_iter->second);
                    if(AvgRefCnvReadMiss != 0 && AvgAltCnvReadMiss != 0){
                        missRateMap[variant.first] = AvgAltCnvReadMiss / (AvgRefCnvReadMiss + AvgAltCnvReadMiss);
                    }
                }
            }
            i++;
        }

    }
}

//filter the variants with high mismatch rate   
void VairiantGraph::filterHighMismatchVariants(std::vector<ReadVariant>& in_readVariant, const Clip& clip, const std::map<int, double>& missRateMap){
    //if the read variant is empty or the cnv is empty or the miss rate map is empty, return
    if(in_readVariant.empty() || clip.cnvVec.empty() || missRateMap.empty()){
        return;
    }

    size_t cnvIndex = 0;

    for(auto& read : in_readVariant){
        //if the read has no variant, continue
        if(read.variantVec.empty()){
            continue;
        }
        auto variantIter = read.variantVec.begin();
        int readStart = read.variantVec.front().position;

        while(cnvIndex > 0 && clip.cnvVec[cnvIndex].first > readStart){
            cnvIndex--;
        }

        //iterate through the variant vector
        while(variantIter != read.variantVec.end()){
            bool shouldErase = false;
            
            //iterate through the cnv vector
            size_t i = cnvIndex;
            while(i < clip.cnvVec.size() && clip.cnvVec[i].first <= variantIter->position){
                //if the variant position is in the range of the cnv, check the miss rate
                if(isPositionInRange(variantIter->position, clip.cnvVec[i].first, clip.cnvVec[i].second)){
                    auto missIter = missRateMap.find(variantIter->position);
                    //if the miss rate is greater than 0.7, erase the variant
                    if(missIter != missRateMap.end() && missIter->second >= 0.7){
                        shouldErase = true;
                        variantIter = read.variantVec.erase(variantIter);
                        break;
                    }
                }
                i++;
            }
            //if the variant should not be erased, increment the variant iterator
            if(!shouldErase){
                ++variantIter;
            }

            cnvIndex = i > 0 ? i - 1 : 0;
        }
    }
}
    
void VairiantGraph::addEdge(std::vector<ReadVariant> &in_readVariant, Clip &clip){

    readVariant = &in_readVariant;
    std::map<std::string,ReadVariant> mergeReadMap;

    // each read will record first and last variant posistion
    std::map<std::string, std::pair<int,int>> alignRange;
    // record an iterator for all alignments of a read.
    std::map<std::string, std::vector<int>> readIdxVec;
    // record need del read index
    std::vector<int> delReadIdx;

    // Check for overlaps among different alignments of a read and filter out the shorter overlapping alignments.
    for (int readIter = 0; readIter < (int)in_readVariant.size(); readIter++) {
        int is_toDelete = 0;
        std::string readName = in_readVariant[readIter].read_name;
        int firstVariantPos = in_readVariant[readIter].variantVec.front().position;
        int lastVariantPos = in_readVariant[readIter].variantVec.back().position;
        auto& readRange = alignRange[readName];
        auto& readIdxVecRef = readIdxVec[readName];

        // Initialize readRange if it's the first appearance
        if (alignRange.find(readName) == alignRange.end()) {
            readRange = {firstVariantPos,lastVariantPos};
        } else {
            // Check for overlaps
            while (readRange.first <= firstVariantPos && firstVariantPos <= readRange.second) {
                if (lastVariantPos < readRange.second) {
                    is_toDelete = 1;
                    delReadIdx.push_back(readIter);
                    break;
                }

                int preAlignIdx = readIdxVecRef.size() - 1;
                if (preAlignIdx < 0 ) break;

                const auto& previousAlignment = in_readVariant[readIdxVecRef[preAlignIdx]];
                const auto& prevVariantVec = previousAlignment.variantVec;
                int prevStart = prevVariantVec.front().position;
                int prevEnd = prevVariantVec.back().position;

                double overlapStart = std::max(prevStart, firstVariantPos);
                double overlapEnd = std::min(prevEnd, lastVariantPos);
                if (overlapStart > overlapEnd) break; // No overlap
                double overlapLen = overlapEnd - overlapStart + 1;

                double alignStart = std::max(prevEnd, lastVariantPos);
                double alignEnd = std::min(prevStart, firstVariantPos);
                double alignSpan = alignStart - alignEnd + 1;
                double overlapRatio = overlapLen / alignSpan;

                // Filtering highly overlapping alignments
                if (overlapRatio >= params->overlapThreshold) {
                    int alignLen1 = prevEnd - prevStart + 1;
                    int alignLen2 = lastVariantPos - firstVariantPos + 1;

                    if (alignLen2 <= alignLen1) {
                        is_toDelete = 1;
                        delReadIdx.push_back(readIter); // Current alignment is shorter
                        break;
                    } else {
                        delReadIdx.push_back(readIdxVecRef[preAlignIdx]); // Previous alignment is shorter
                        readIdxVecRef.pop_back();
                        readRange.second = (preAlignIdx > 0) ? in_readVariant[readIdxVecRef[preAlignIdx - 1]].variantVec.back().position : firstVariantPos;
                    }
                } else {
                    break;
                }
            }
            // update range
            readRange.second = lastVariantPos;
        }
        if (is_toDelete == 0 )
            readIdxVecRef.push_back(readIter);
    }

    // sort read index
    std::sort(delReadIdx.begin(), delReadIdx.end());
    // remove overlap alignment
    delReadIdx.push_back((int)in_readVariant.size());
    int saveIter = *(delReadIdx.begin());
    for (auto delIter = delReadIdx.begin(), nextdelIter = std::next(delReadIdx.begin(), 1); nextdelIter != delReadIdx.end(); delIter++ , nextdelIter++) {
        auto nowDelIter = *delIter+1;
        while (nowDelIter<*nextdelIter){
            in_readVariant[saveIter++]=in_readVariant[nowDelIter++];
        }
    }
    in_readVariant.erase( std::next(in_readVariant.begin(), saveIter), in_readVariant.end());
    
    CnvStatistics cnvStats;
    //calculate the mismatch rate of the cnv
    calculateCnvMismatchRate(in_readVariant, clip);
    //aggregate the mismatch rate of the cnv
    aggregateCnvReadMismatchRate(in_readVariant, clip, cnvStats.cnvReadMmrate);
    //calculate the average mismatch rate of the cnv
    calculateAverageMismatchRate(clip, cnvStats.cnvReadMmrate, cnvStats.missRateMap);
    //filter the variants with high mismatch rate
    filterHighMismatchVariants(in_readVariant, clip, cnvStats.missRateMap);

    int readCount=0;
    // merge alignment
    for(std::vector<ReadVariant>::iterator readIter = in_readVariant.begin() ; readIter != in_readVariant.end() ; readIter++ ){
 
        // Creating a pseudo read which allows filtering out variants that should not be phased
        //ReadVariant tmpRead;
        // Visiting all the variants on the read
        for( auto variant : (*readIter).variantVec ){
            readCount++;
            // modification
            if( variant.quality == -2 || variant.quality == -3 ){
                (*variantType)[variant.position] = 2;
                variant.quality = 60;
            }
            // structure variation
            else if( variant.quality == -1 ){
                (*variantType)[variant.position] = 1;
                if( variant.allele == 1 ){
                    // SVcaller calling
                    variant.quality = 60; 
                }
                else{
                    // In SVcaller, unmarked reads are assumed to be REF
                    variant.quality = 30;
                }
            }
            // indel
            else if( variant.quality == -4 ){
                (*variantType)[variant.position] = 3;
                variant.quality = 60;
            }
            //danger indel
            else if( variant.quality == -5 ){
                (*variantType)[variant.position] = 4;
                variant.quality = 60;
            }
            // The remaining variants will be labeled as SNPs
            else{
                (*variantType)[variant.position] = 0;
            }
            mergeReadMap[(*readIter).read_name].variantVec.push_back(variant);

            //tmpRead.variantVec.push_back(variant);
            
            // Each position will record the included reads and their corresponding base qualities.
            auto variantIter = totalVariantInfo->find(variant.position);
            
            if( variantIter == totalVariantInfo->end() ){
                (*totalVariantInfo)[variant.position] = new ReadBaseMap();
            }
            
            (*(*totalVariantInfo)[variant.position])[(*readIter).read_name] = variant.quality;
        }
    }   

    for(std::map<std::string,ReadVariant>::iterator readIter = mergeReadMap.begin() ; readIter != mergeReadMap.end() ; readIter++){ 
        (*readIter).second.sort();
        
        // iter all pair of snp and construct initial graph
        std::vector<Variant>::iterator variant1Iter = readIter->second.variantVec.begin();
        std::vector<Variant>::iterator variant2Iter = std::next(variant1Iter,1);
        
        while(variant1Iter != readIter->second.variantVec.end() && variant2Iter != readIter->second.variantVec.end() ){
            // create new edge if not exist
            std::map<int,VariantEdge*>::iterator posIter = edgeList->find(variant1Iter->position);
            if( posIter == edgeList->end() )
                (*edgeList)[variant1Iter->position] = new VariantEdge(variant1Iter->position);

            // add edge process
            for(int nextNode = 0 ; nextNode < params->connectAdjacent; nextNode++){
                // this allele support ref
                if( variant1Iter->allele == 0 )
                    (*edgeList)[variant1Iter->position]->ref->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).first,params->baseQuality,params->edgeWeight);
                // this allele support alt
                if( (*variant1Iter).allele == 1 )
                    (*edgeList)[variant1Iter->position]->alt->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).first,params->baseQuality,params->edgeWeight);
                
                // next snp
                variant2Iter++;
                if( variant2Iter == readIter->second.variantVec.end() ){
                    break;
                }
            }

            variant1Iter++;
            variant2Iter = std::next(variant1Iter,1);
        }

        //count the ref and alt base amount of the last variant on the read
	    if ( variant1Iter != (*readIter).second.variantVec.end() && variant2Iter == (*readIter).second.variantVec.end() ) {
            std::map<int,VariantEdge*>::iterator posIter = edgeList->find((*variant1Iter).position);
            if( posIter == edgeList->end() ) {
                (*edgeList)[(*variant1Iter).position] = new VariantEdge((*variant1Iter).position);
            }
        }
    }
} 

void VairiantGraph::readCorrection(){

    
    std::map<std::string,std::map<int,std::map<int,int>>> readBlockHP;
    
    std::map<std::string,std::map<int,std::map<int,int>>> readBlockHPcount;
    
    
    // haplotype, <position <allele, base count>>
    std::map<int,std::map<int,std::map<double,double>>> *hpAlleleCountMap = new std::map<int,std::map<int,std::map<double,double>>>;

    
    // iter all read, determine the haplotype of the read
    for(std::vector<ReadVariant>::iterator readIter = (*readVariant).begin() ; readIter != (*readVariant).end() ; readIter++ ){
        double refCount = 0;
        double altCount = 0;
        //int block;  
          
        // loop all variant 
        for( auto variant : (*readIter).variantVec ){
            PosAllele refAllele = std::make_pair( variant.position , variant.allele+1);
            std::map<PosAllele,int>::iterator nodePS = bkResult->find(refAllele);
        
            //block = nodePS->second;
            if( nodePS != bkResult->end() ){
                if((*bkResult)[refAllele] != 0 ){
                    if((*variantType)[variant.position] == 0){
                        if((*subNodeHP)[refAllele]==0)refCount++;
                        else altCount++;
                    }
                    else if((*variantType)[variant.position] == 1){
                        if((*subNodeHP)[refAllele]==0)refCount++;
                        else altCount++;
                    }
                    else if((*variantType)[variant.position] == 2){
                        continue;
                    }
                    else if((*variantType)[variant.position] == 3){
                        if((*subNodeHP)[refAllele]==0)refCount+=0.1;
                        else altCount+=0.1;
                    }
                    else if((*variantType)[variant.position] == 4){
                        if((*subNodeHP)[refAllele]==0)refCount+=0.1;
                        else altCount+=0.1;
                    }
                }
            }

        }

        // tag high confident reads
        if( std::max(refCount,altCount)/(refCount+altCount) > params->readConfidence && (refCount + altCount) > 1 ){
            // tag read with the corresponding haplotype
            int belongHP = ( refCount > altCount ? 0 : 1 );
            (*readHpMap)[(*readIter).read_name] = belongHP;
            
            //readBlockHP[(*readIter).read_name][(*readIter).reference_start][block]=belongHP;
            //readBlockHPcount[(*readIter).read_name][block][belongHP]++;

            for(auto variantIter = (*readIter).variantVec.begin() ; variantIter != (*readIter).variantVec.end() ; variantIter++ ){
                if( (*variantIter).allele == 0 || (*variantIter).allele == 1){
                    (*hpAlleleCountMap)[belongHP][(*variantIter).position][(*variantIter).allele]++;
                }
            }
        }
        else{
            (*readHpMap)[(*readIter).read_name] = -1;
        }
    }

    /*
    for(auto readIter = readBlockHP.begin() ; readIter != readBlockHP.end() ; readIter++ ){
        
        int max = 0;
        for(auto blockIter = readBlockHPcount[readIter->first].begin() ; blockIter != readBlockHPcount[readIter->first].end() ; blockIter++ ){
            if( (*blockIter).second.size() > max ){
                max = (*blockIter).second.size();
            }
        }

        std::cout<< readIter->first << "\t" << max;
        
        for(auto startIter = readIter->second.begin() ; startIter != readIter->second.end() ; startIter++ ){
            std::cout<< "\t" << startIter->first << ",";
            for(auto blockIter = startIter->second.begin() ; blockIter != startIter->second.end() ; blockIter++ ){
                std::cout<< blockIter->first << "," <<  blockIter->second;
            }
        }
        std::cout<< "\n";
    }
    */

    double snpConfidenceThreshold = params->snpConfidence;

    subNodeHP->clear();
    
    std::map<int,std::map<int,int>> hpAllele;
    // reassign allele result
    for(auto variantIter = totalVariantInfo->begin() ; variantIter != totalVariantInfo->end() ; variantIter++ ){
        int position = variantIter->first;
        PosAllele refAllele = std::make_pair(position, 1);
        PosAllele altAllele = std::make_pair(position, 2);
        
        double hp1Ref = (*hpAlleleCountMap)[0][position][0];
        double hp1Alt = (*hpAlleleCountMap)[0][position][1];
        double hp2Ref = (*hpAlleleCountMap)[1][position][0];
        double hp2Alt = (*hpAlleleCountMap)[1][position][1];
        double result1reads = hp1Ref + hp2Alt;
        double result2reads = hp2Ref + hp1Alt;
        double resultConfidence = std::max(result1reads, result2reads) / (result1reads + result2reads);
        
        int hp1Result = -1;
        int hp2Result = -1;
        
        //std::cout << "RC\t" << position+1 << "\t" << result1reads << "\t" << result2reads << "\t" << resultConfidence << "\n" ;
        
        if( resultConfidence > snpConfidenceThreshold ){
            if( result1reads > result2reads ){
                hp1Result = 0;
                hp2Result = 1;
            }
            else if( result1reads < result2reads ){
                hp1Result = 1;
                hp2Result = 0;
            }
        }

        if( hp1Result != -1 && hp2Result != -1 ){
            (*subNodeHP)[refAllele] = hp1Result;
            (*subNodeHP)[altAllele] = hp2Result;
        }
        else{
            bkResult->erase(refAllele);
            bkResult->erase(altAllele);
        }
    }

    delete hpAlleleCountMap;
}

void VairiantGraph::writingDotFile(std::string dotPrefix){
    
    std::ofstream resultVcf(dotPrefix+".dot");

    if(!resultVcf.is_open()){
        std::cerr<< "Fail to open write file: " << dotPrefix+".vcf" << "\n";
    }
    else{
        resultVcf << "digraph G {\n";

        for(auto edge : dotResult){
            resultVcf << edge << "\n";
        }
        resultVcf << "}\n";
    }
    return;
}

void VairiantGraph::exportResult(std::string chrName, PhasingResult &result){
    
    // loop all position
    for( std::map<int,ReadBaseMap*>::iterator variantIter = totalVariantInfo->begin() ; variantIter != totalVariantInfo->end() ; variantIter++ ){
        
        PhasingElement tmp;
        
        PosAllele ref = std::make_pair( variantIter->first , 1);
        PosAllele alt = std::make_pair( variantIter->first , 2);
        
        std::map<PosAllele,int>::iterator psRefIter = bkResult->find(ref);
        std::map<PosAllele,int>::iterator psAltIter = bkResult->find(alt);
        
        if( psRefIter != bkResult->end() || psAltIter != bkResult->end() ){
            if( psRefIter != bkResult->end() )
                tmp.block = (*psRefIter).second;
            else
                tmp.block = (*psAltIter).second;
            tmp.RAstatus = std::to_string((*subNodeHP)[ref]) + "|" + std::to_string((*subNodeHP)[alt]);
        }
        else
            continue;
        
        if( tmp.block != 0){
            std::string key = chrName + "_" + std::to_string( variantIter->first );
            result[key] = tmp;
        }
    }
}

std::map<std::string,int>* VairiantGraph::getReadHP(){
    return readHpMap;
}

int VairiantGraph::totalNode(){
    return totalVariantInfo->size();
}

void VairiantGraph::phasingProcess(){
    // This step involves converting all reads into a graph structure, which will be stored as an edge list
    // in a two-layer map. The first layer of the map uses the starting coordinate as the key and contains
    // a second layer map as the value. The second layer map uses the destination coordinate as the key and
    // stores the number of support read as values. (There is another map used for debugging purposes that
    // treats the read name vector as a value.) The method begins by visiting the coordinates covered by each 
    // read and recording this information in 'totalVariantInfo.' Subsequently, it connects the coordinates contained 
    // in each read on the graph. Specifically, each coordinate is connected to the next N coordinates in a 
    // linear fashion.
    this->edgeConnectResult();

    // This step will utilize the results of graph phasing to attempt to separate all the reads into two 
    // haplotypes and then identify high-confidence SNPs using reads from the two distinct haplotypes.
    this->readCorrection();  
}

Clip::Clip(std::string &chr, ClipCount &clipCount){
    this->chr = chr;
    getCNVInterval(clipCount);
}

Clip::~Clip(){
}

//update the threshold
void Clip::updateThreshold(int upCount){
    state.rejectCount = upCount;
    if(upCount >= 20){
        state.pullDownCount = upCount / 2;
        state.slowDownCount = 5;
    }
    else if(upCount >= 10){
        state.pullDownCount = upCount / 2;
        state.slowDownCount = upCount / 4;
    }
    else{
        state.pullDownCount = 5;
        state.slowDownCount = 2;
    }
}

void Clip::getCNVInterval(ClipCount &clipCount){
    int upCount = 0;
    int downCount = 0;
    int AreaSize = 30000;
    state.reset();

    clipCount[clipCount.rbegin()->first + AreaSize] = clipCount.rbegin()->second;
    std::map<std::string, std::map<int,int>> cnvArea;


    for(auto posIter = clipCount.begin(); posIter != clipCount.end() ; posIter++ ){
        upCount = posIter->second[FRONT];
        downCount = posIter->second[BACK];

        //if the current state is not push, not slowdown, and not slowup
        if(!state.push && !state.slowDown && !state.slowUp){
            //if the up count is greater than 5 and the current count is 0, then change the state to push and slow down
            if(upCount >= 5 && state.currCount == 0){
                state.push = 1;
                state.slowUp = 0;
                state.slowDown = 1;
                state.currCount = upCount - downCount;
                state.candidateStartPos = posIter->first;
                state.candidateEndPos = posIter->first + AreaSize;
                updateThreshold(upCount);
            }
            //if the up count is greater than the down count and the current count is 0, then change the state to slowup
            else if(upCount > downCount && state.currCount == 0){
                state.push = 0;
                state.slowUp = 1;
                state.slowDown = 0;
                state.currCount = upCount - downCount;
                state.candidateStartPos = posIter->first;
                state.candidateEndPos = posIter->first + AreaSize;
            }
        }
        //if the current state is push and slowdown
        else if(state.push && state.slowDown){
            //if the up count is greater than the reject count, then change the state to push
            if(upCount > state.rejectCount){
                state.push = 1;
                state.slowUp = 0;
                state.slowDown = 1;
                updateThreshold(upCount);
                state.candidateStartPos = posIter->first;
                state.candidateEndPos = posIter->first + AreaSize;
            }
            //update the current count
            state.currCount = state.currCount + upCount - downCount;
            //if the current count is greater than 30, then change the candidate end position to the current position + AreaSize
            if(state.currCount > 30){
                state.candidateEndPos = posIter->first + AreaSize;
            }
            //if the down count is greater than the pull down count, then push the cnv to the vector
            if(downCount >= state.pullDownCount){
                cnvVec.emplace_back(state.candidateStartPos, posIter->first);
                state.reset();
            }
            //if the current count is less than or equal to the slow down count and the current position is less than or equal to the candidate end position, then push the cnv to the vector
            else if(state.currCount <= state.slowDownCount && posIter->first <= state.candidateEndPos){
                cnvVec.emplace_back(state.candidateStartPos, posIter->first);            
                state.reset();
            }
            //if the current position is greater than the candidate end position or the current count is less than or equal to 0 or the current position is greater than or equal to the candidate start position + 200000, then reset the state
            if(posIter->first > state.candidateEndPos || state.currCount <= 0 || posIter->first - state.candidateStartPos >= 200000){
                state.reset();
            }
        }
        //if the current state is slow up
        else if(state.slowUp){
            //if the current count is greater than 20, then down count is greater than the current count / 4, then push the cnv to the vector
            if(state.currCount > 20 ? downCount >= state.currCount/4 : downCount >= 5){
                cnvVec.emplace_back(state.candidateStartPos, posIter->first);
                state.reset();
            }
            //if the up count is greater than 5, then change the state to push and slow down
            else if(upCount >= 5){
                state.push = 1;
                state.slowUp = 0;
                state.slowDown = 1;
                state.currCount = upCount - downCount;
                state.candidateStartPos = posIter->first;
                state.candidateEndPos = posIter->first + AreaSize;
                updateThreshold(upCount);
            }
            else{
                state.currCount = state.currCount + upCount - downCount;
                //if the current count is greater than 30, then change the candidate end position to the current position + AreaSize
                if(state.currCount > 30){
                    state.candidateEndPos = posIter->first + AreaSize;
                }
                //if the current position is greater than the candidate end position or the current count is less than or equal to 0 or the current position is greater than or equal to the candidate start position + 200000, then reset the state
                if(posIter->first > state.candidateEndPos || state.currCount <= 0 || posIter->first - state.candidateStartPos >= 200000){
                    state.reset();
                }
            }
        }
    }
    clipCount.erase(--clipCount.end());
}