#include "PhasingGraph.h"

namespace {

void appendDotEdges(std::vector<std::string>& outDotResult, int currPos, const std::pair<PosAllele, PosAllele>& edgePair){
    std::string refEdge = std::to_string(currPos + 1) + ".1\t->\t" + std::to_string(edgePair.first.first + 1) + "." + std::to_string(edgePair.first.second);
    std::string altEdge = std::to_string(currPos + 1) + ".2\t->\t" + std::to_string(edgePair.second.first + 1) + "." + std::to_string(edgePair.second.second);

    outDotResult.push_back(refEdge);
    outDotResult.push_back(altEdge);
}

void recordVoteForNextPosition(
    int currHP,
    int targetPos,
    int targetHaplotypeRelation,
    VoteResult& vote,
    std::map<int, std::map<int, float> >& hpCountMap2,
    std::map<int, std::vector<VoteResult> >& hpCountMap3){
    if(currHP == 1){
        if(targetHaplotypeRelation == 1){
            hpCountMap2[targetPos][1] += vote.weight;
            vote.hap = 1;
        }
        else if(targetHaplotypeRelation == 2){
            hpCountMap2[targetPos][2] += vote.weight;
            vote.hap = 2;
        }
    }
    else if(currHP == 2){
        if(targetHaplotypeRelation == 1){
            hpCountMap2[targetPos][2] += vote.weight;
            vote.hap = 2;
        }
        else if(targetHaplotypeRelation == 2){
            hpCountMap2[targetPos][1] += vote.weight;
            vote.hap = 1;
        }
    }

    hpCountMap3[targetPos].push_back(vote);
}

}
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
void VairiantGraph::scanVariantsAndBuildBlocks(std::map<int, int>& hpResult, PhasedBlocks& phasedBlocks, std::vector<std::string>& outDotResult){
    std::map<int, std::vector<VoteResult> > hpCountMap3;
    std::map<int, std::map<int,float> > hpCountMap2;
    int blockStart = -1;
    int lastConnectPos = -1;

    for(std::map<int,ReadBaseMap*>::iterator variantIter = totalVariantInfo->begin(); variantIter != totalVariantInfo->end(); variantIter++ ){
        std::map<int,ReadBaseMap*>::iterator nextNodeIter = std::next(variantIter, 1);
        if(nextNodeIter == totalVariantInfo->end()){
            break;
        }

        int currPos = variantIter->first;
        int nextPos = nextNodeIter->first;

        if(std::abs(nextPos - currPos) > params->distance){
            continue;
        }

        float h1 = hpCountMap2[currPos][1];
        float h2 = hpCountMap2[currPos][2];

        std::pair<float, float> special = Onelongcase(hpCountMap3[currPos]);
        if(special.first != -1){
            h1 = special.first;
            h2 = special.second;
        }

        if(h1 == h2){
            if(currPos < lastConnectPos){
                continue;
            }

            blockStart = currPos;
            phasedBlocks[blockStart].push_back(currPos);
            hpResult[currPos] = 1;
        }
        else{
            int currHP = (h1 > h2 ? 1 : 2);
            hpResult[currPos] = currHP;
            phasedBlocks[blockStart].push_back(currPos);
        }

        std::map<int,VariantEdge*>::iterator edgeIter = edgeList->find(currPos);
        if(edgeIter == edgeList->end()){
            continue;
        }

        for(int i = 0; i < params->connectAdjacent; i++ ){
            VoteResult vote(currPos, 1);
            std::pair<PosAllele,PosAllele> bestEdgePair = edgeIter->second->findBestEdgePair(nextNodeIter->first, params->isONT, params->edgeThreshold, false, *variantType, vote);

            if((*variantType)[currPos] == 4){
                vote.weight = 0.1;
            }

            if(bestEdgePair.first.second != -1){
                recordVoteForNextPosition(
                    hpResult[currPos],
                    nextNodeIter->first,
                    bestEdgePair.first.second,
                    vote,
                    hpCountMap2,
                    hpCountMap3);

                if(params->generateDot){
                    appendDotEdges(outDotResult, currPos, bestEdgePair);
                }

                lastConnectPos = nextNodeIter->first;
            }

            nextNodeIter++;
            if(nextNodeIter == totalVariantInfo->end()){
                break;
            }
        }
    }
}

void VairiantGraph::materializeBlockResults(
    const std::map<int, int>& hpResult,
    const PhasedBlocks& phasedBlocks,
    std::map<PosAllele, int>& outBkResult,
    std::map<PosAllele, int>& outSubNodeHP){
    for(auto blockIter = phasedBlocks.begin(); blockIter != phasedBlocks.end(); blockIter++ ){
        if((*blockIter).second.size() <= 1){
            continue;
        }

        for(auto currIter = (*blockIter).second.begin(); currIter != (*blockIter).second.end(); currIter++ ){
            auto nextIter = std::next(currIter, 1);
            if(nextIter == (*blockIter).second.end()){
                continue;
            }

            PosAllele refStart = std::make_pair((*currIter), 1);
            PosAllele altStart = std::make_pair((*currIter), 2);
            PosAllele refEnd = std::make_pair((*nextIter), 1);
            PosAllele altEnd = std::make_pair((*nextIter), 2);

            outBkResult[refStart] = (*blockIter).first + 1;
            outBkResult[refEnd]   = (*blockIter).first + 1;
            outBkResult[altStart] = (*blockIter).first + 1;
            outBkResult[altEnd]   = (*blockIter).first + 1;

            if(currIter == (*blockIter).second.begin()){
                outSubNodeHP[refStart] = 0;
                outSubNodeHP[altStart] = 1;
            }

            if(hpResult.at(*currIter) == 0 || hpResult.at(*nextIter) == 0){

            }
            else if(hpResult.at(*currIter) == hpResult.at(*nextIter)){
                outSubNodeHP[refEnd] = outSubNodeHP[refStart];
                outSubNodeHP[altEnd] = outSubNodeHP[altStart];
            }
            else{
                outSubNodeHP[refEnd] = outSubNodeHP[altStart];
                outSubNodeHP[altEnd] = outSubNodeHP[refStart];
            }
        }
    }
}

void VairiantGraph::edgeConnectResult(){
    std::map<int, int> hpResult;
    PhasedBlocks phasedBlocks;

    scanVariantsAndBuildBlocks(hpResult, phasedBlocks, dotResult);
    materializeBlockResults(hpResult, phasedBlocks, *bkResult, *subNodeHP);
}

VairiantGraph::VairiantGraph(std::string &in_ref, PhasingParameters &in_params, std::string &in_chrName){
    params=&in_params;
    ref=&in_ref;
    totalVariantInfo = new std::map<int,ReadBaseMap*>;
    edgeList = new std::map<int,VariantEdge*>;
    bkResult = new std::map<PosAllele,int>;
    subNodeHP = new std::map<PosAllele,int>;
    variantType = new std::map<int,int>;
    readHpMap = new std::map<std::string,int>;
    chrName = &in_chrName;
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
    
void VairiantGraph::filterOverlappingAlignments(std::vector<ReadVariant> &in_readVariant){
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

        // On the first appearance, operator[] initializes readRange to {0, 0},
        // so the overlap condition below naturally does not hold.
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
}

void VairiantGraph::applyCnvFilter(std::vector<ReadVariant> &in_readVariant, Clip &clip){
    CnvStatistics cnvStats;
    //calculate the mismatch rate of the cnv
    calculateCnvMismatchRate(in_readVariant, clip);
    //aggregate the mismatch rate of the cnv
    aggregateCnvReadMismatchRate(in_readVariant, clip, cnvStats.cnvReadMmrate);
    //calculate the average mismatch rate of the cnv
    calculateAverageMismatchRate(clip, cnvStats.cnvReadMmrate, cnvStats.missRateMap);
    //filter the variants with high mismatch rate
    filterHighMismatchVariants(in_readVariant, clip, cnvStats.missRateMap);
}

void VairiantGraph::buildVariantGraph(std::vector<ReadVariant> &in_readVariant){
    std::map<std::string,ReadVariant> mergeReadMap;

    // merge alignment
    for(std::vector<ReadVariant>::iterator readIter = in_readVariant.begin() ; readIter != in_readVariant.end() ; readIter++ ){

        // Visiting all the variants on the read
        for( auto variant : (*readIter).variantVec ){
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

void VairiantGraph::addEdge(std::vector<ReadVariant> &in_readVariant, Clip &clip){

    readVariant = &in_readVariant;

    filterOverlappingAlignments(in_readVariant);
    applyCnvFilter(in_readVariant, clip);
    buildVariantGraph(in_readVariant);

} 

void VairiantGraph::accumulateReadHaplotypeEvidence(const Variant& variant, double& refCount, double& altCount){
    PosAllele refAllele = std::make_pair(variant.position, variant.allele + 1);
    std::map<PosAllele, int>::iterator nodePS = bkResult->find(refAllele);

    if(nodePS == bkResult->end() || nodePS->second == 0){
        return;
    }

    int refHaplotype = (*subNodeHP)[refAllele];
    int positionVariantType = (*variantType)[variant.position];

    if(positionVariantType == 0 || positionVariantType == 1){
        if(refHaplotype == 0) refCount++;
        else altCount++;
    }
    else if(positionVariantType == 2){
        return;
    }
    else if(positionVariantType == 3 || positionVariantType == 4){
        if(refHaplotype == 0) refCount += 0.1;
        else altCount += 0.1;
    }
}

void VairiantGraph::assignReadHaplotypesAndCollectAlleleCounts(HpAlleleCountMap& hpAlleleCountMap){
    // iter all read, determine the haplotype of the read
    for(std::vector<ReadVariant>::iterator readIter = (*readVariant).begin(); readIter != (*readVariant).end(); readIter++){
        double refCount = 0;
        double altCount = 0;

        // loop all variant
        for(const auto& variant : (*readIter).variantVec){
            accumulateReadHaplotypeEvidence(variant, refCount, altCount);
        }

        // tag high confident reads
        if(std::max(refCount, altCount) / (refCount + altCount) > params->readConfidence && (refCount + altCount) > 1){
            // tag read with the corresponding haplotype
            int belongHP = (refCount > altCount ? 0 : 1);
            (*readHpMap)[(*readIter).read_name] = belongHP;

            for(auto variantIter = (*readIter).variantVec.begin(); variantIter != (*readIter).variantVec.end(); variantIter++){
                if((*variantIter).allele == 0 || (*variantIter).allele == 1){
                    hpAlleleCountMap[belongHP][(*variantIter).position][(*variantIter).allele]++;
                }
            }
        }
        else{
            (*readHpMap)[(*readIter).read_name] = -1;
        }
    }
}

void VairiantGraph::reassignVariantHaplotypes(const HpAlleleCountMap& hpAlleleCountMap){
    auto getHpAlleleCount = [&](int haplotype, int position, int allele){
        auto hpIter = hpAlleleCountMap.find(haplotype);
        if(hpIter == hpAlleleCountMap.end()){
            return 0.0;
        }

        auto positionIter = hpIter->second.find(position);
        if(positionIter == hpIter->second.end()){
            return 0.0;
        }

        auto alleleIter = positionIter->second.find(allele);
        if(alleleIter == positionIter->second.end()){
            return 0.0;
        }

        return alleleIter->second;
    };

    double snpConfidenceThreshold = params->snpConfidence;

    subNodeHP->clear();

    // reassign allele result
    for(auto variantIter = totalVariantInfo->begin(); variantIter != totalVariantInfo->end(); variantIter++){
        int position = variantIter->first;
        PosAllele refAllele = std::make_pair(position, 1);
        PosAllele altAllele = std::make_pair(position, 2);

        double hp1Ref = getHpAlleleCount(0, position, 0);
        double hp1Alt = getHpAlleleCount(0, position, 1);
        double hp2Ref = getHpAlleleCount(1, position, 0);
        double hp2Alt = getHpAlleleCount(1, position, 1);
        double result1reads = hp1Ref + hp2Alt;
        double result2reads = hp2Ref + hp1Alt;
        double resultConfidence = std::max(result1reads, result2reads) / (result1reads + result2reads);

        int hp1Result = -1;
        int hp2Result = -1;

        if(resultConfidence > snpConfidenceThreshold){
            if(result1reads > result2reads){
                hp1Result = 0;
                hp2Result = 1;
            }
            else if(result1reads < result2reads){
                hp1Result = 1;
                hp2Result = 0;
            }
        }

        if(hp1Result != -1 && hp2Result != -1){
            (*subNodeHP)[refAllele] = hp1Result;
            (*subNodeHP)[altAllele] = hp2Result;
        }
        else{
            bkResult->erase(refAllele);
            bkResult->erase(altAllele);
        }
    }
}

void VairiantGraph::readCorrection(){
    HpAlleleCountMap hpAlleleCountMap;
    assignReadHaplotypesAndCollectAlleleCounts(hpAlleleCountMap);
    reassignVariantHaplotypes(hpAlleleCountMap);
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

void Clip::handleIdleState(int upCount, int downCount, int pos){
    if(upCount >= CNV_PUSH_THRESHOLD && state.currCount == 0){
        state.push = 1;
        state.slowUp = 0;
        state.slowDown = 1;
        state.currCount = upCount - downCount;
        state.candidateStartPos = pos;
        state.candidateEndPos = pos + CNV_AREA_SIZE;
        updateThreshold(upCount);
    }
    else if(upCount > downCount && state.currCount == 0){
        state.push = 0;
        state.slowUp = 1;
        state.slowDown = 0;
        state.currCount = upCount - downCount;
        state.candidateStartPos = pos;
        state.candidateEndPos = pos + CNV_AREA_SIZE;
    }
}

void Clip::handleActiveState(int upCount, int downCount, int pos){
    if(upCount > state.rejectCount){
        state.push = 1;
        state.slowUp = 0;
        state.slowDown = 1;
        updateThreshold(upCount);
        state.candidateStartPos = pos;
        state.candidateEndPos = pos + CNV_AREA_SIZE;
    }
    state.currCount = state.currCount + upCount - downCount;
    if(state.currCount > CNV_EXTEND_THRESHOLD){
        state.candidateEndPos = pos + CNV_AREA_SIZE;
    }
    if(downCount >= state.pullDownCount){
        cnvVec.emplace_back(state.candidateStartPos, pos);
        state.reset();
    }
    else if(state.currCount <= state.slowDownCount && pos <= state.candidateEndPos){
        cnvVec.emplace_back(state.candidateStartPos, pos);
        state.reset();
    }
    if(pos > state.candidateEndPos || state.currCount <= 0 || pos - state.candidateStartPos >= CNV_MAX_LENGTH){
        state.reset();
    }
}

void Clip::handleSlowUpState(int upCount, int downCount, int pos){
    if(state.currCount > CNV_SLOWUP_COMMIT_THRESHOLD ? downCount >= state.currCount/4 : downCount >= CNV_PUSH_THRESHOLD){
        cnvVec.emplace_back(state.candidateStartPos, pos);
        state.reset();
    }
    else if(upCount >= CNV_PUSH_THRESHOLD){
        state.push = 1;
        state.slowUp = 0;
        state.slowDown = 1;
        state.currCount = upCount - downCount;
        state.candidateStartPos = pos;
        state.candidateEndPos = pos + CNV_AREA_SIZE;
        updateThreshold(upCount);
    }
    else{
        state.currCount = state.currCount + upCount - downCount;
        if(state.currCount > CNV_EXTEND_THRESHOLD){
            state.candidateEndPos = pos + CNV_AREA_SIZE;
        }
        if(pos > state.candidateEndPos || state.currCount <= 0 || pos - state.candidateStartPos >= CNV_MAX_LENGTH){
            state.reset();
        }
    }
}

void Clip::getCNVInterval(ClipCount &clipCount){
    if (clipCount.empty()) {
        return;
    }
    state.reset();

    clipCount[clipCount.rbegin()->first + CNV_AREA_SIZE] = clipCount.rbegin()->second;

    for(auto posIter = clipCount.begin(); posIter != clipCount.end(); posIter++){
        int upCount   = posIter->second[FRONT];
        int downCount = posIter->second[BACK];
        int pos       = posIter->first;

        if(!state.push && !state.slowDown && !state.slowUp)
            handleIdleState(upCount, downCount, pos);
        else if(state.push && state.slowDown)
            handleActiveState(upCount, downCount, pos);
        else if(state.slowUp)
            handleSlowUpState(upCount, downCount, pos);
    }
    clipCount.erase(--clipCount.end());
}
