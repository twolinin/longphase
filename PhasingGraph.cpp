#include "PhasingGraph.h"

#include <cmath>   // std::log2 for Shannon-entropy computation
#include <cstdlib> // std::strtof

namespace {

void appendDotEdges(std::string& outDotResult, int currPos, const std::pair<PosAllele, PosAllele>& edgePair, float h1Weight, float h2Weight){
    // The label on each edge is the accumulated vote weight seen at
    // currPos on that haplotype (HP1 for the .1 edge, HP2 for .2).
    // These labels are new in v2.1 (earlier versions wrote unlabelled
    // edges). The GNN gets the same weights via appendGnnEdges.
    std::string refEdge = std::to_string(currPos + 1) + ".1\t->\t" + std::to_string(edgePair.first.first + 1) + "." + std::to_string(edgePair.first.second) + "\t[label=" + std::to_string(h1Weight) + "]";
    std::string altEdge = std::to_string(currPos + 1) + ".2\t->\t" + std::to_string(edgePair.second.first + 1) + "." + std::to_string(edgePair.second.second) + "\t[label=" + std::to_string(h2Weight) + "]";

    outDotResult += refEdge;
    outDotResult += '\n';
    outDotResult += altEdge;
    outDotResult += '\n';
}

// The weight of the GNN edges. The model was trained on weights read back
// from the DOT label, so they take the same text round trip
// (std::to_string, then parsed as float) to give identical values.
float gnnEdgeWeight(float weight){
    return std::strtof(std::to_string(weight).c_str(), NULL);
}

void recordVoteForNextPosition(
    int currHP,
    int targetPos,   // variant index
    int targetHaplotypeRelation,
    VoteResult& vote,
    std::vector<std::array<float, 2> >& hpCountMap2,
    std::vector<std::vector<VoteResult> >& hpCountMap3){
    if(currHP == 1){
        if(targetHaplotypeRelation == 1){
            hpCountMap2[targetPos][0] += vote.weight;
            vote.hap = 1;
        }
        else if(targetHaplotypeRelation == 2){
            hpCountMap2[targetPos][1] += vote.weight;
            vote.hap = 2;
        }
    }
    else if(currHP == 2){
        if(targetHaplotypeRelation == 1){
            hpCountMap2[targetPos][1] += vote.weight;
            vote.hap = 2;
        }
        else if(targetHaplotypeRelation == 2){
            hpCountMap2[targetPos][0] += vote.weight;
            vote.hap = 1;
        }
    }

    hpCountMap3[targetPos].push_back(vote);
}

}
//SubEdge

SubEdge::SubEdge():readCount(0), refHint(0), altHint(0){ 
}

SubEdge::~SubEdge(){ 
}

void SubEdge::destroy(){
    CountTable().swap(refReadCount);
    CountTable().swap(altReadCount);
}

void SubEdge::shrinkToFit(){
    refReadCount.shrink_to_fit();
    altReadCount.shrink_to_fit();
}

// Count of pos in table, inserted as 0 if absent (as std::map::operator[]). The positions of one
// read arrive in increasing order, so the search starts at the last position touched.
float &SubEdge::findOrInsert(CountTable &table, size_t &hint, int pos){
    auto less = [](const std::pair<int, float> &entry, int key){ return entry.first < key; };
    CountTable::iterator it;
    if( hint < table.size() && table[hint].first <= pos )
        it = std::lower_bound(table.begin() + hint, table.end(), pos, less);
    else
        it = std::lower_bound(table.begin(), table.begin() + std::min(hint, table.size()), pos, less);
    if( it == table.end() || it->first != pos )
        it = table.insert(it, std::make_pair(pos, 0.0f));
    hint = it - table.begin();
    return it->second;
}

float SubEdge::findCount(const CountTable &table, int pos){
    auto less = [](const std::pair<int, float> &entry, int key){ return entry.first < key; };
    CountTable::const_iterator it = std::lower_bound(table.begin(), table.end(), pos, less);
    if( it != table.end() && it->first == pos )
        return it->second;
    return 0;
}

void SubEdge::addSubEdge(int currentQuality, const Variant &connectNode, const std::string &readName, int baseQuality, double edgeWeight){
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
        float &count = findOrInsert(refReadCount, refHint, connectNode.position);
        if ( currentQuality >= baseQuality && connectNode.quality >= baseQuality)
            count++;
        else
            count = count + edgeWeight ;
        
        
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
        float &count = findOrInsert(altReadCount, altHint, connectNode.position);
        if ( currentQuality >= baseQuality && connectNode.quality >= baseQuality)
            count++;
        else
            count = count + edgeWeight ;
        
	//(*altReadCount)[connectNode.position]++;
    }
    readCount++;
}

std::pair<float,float> SubEdge::BestPair(int targetPos){
    return std::make_pair( getRefReadCount(targetPos), getAltReadCount(targetPos) );
}

float SubEdge::getRefReadCount(int targetPos){
    return findCount(refReadCount, targetPos);
}

float SubEdge::getAltReadCount(int targetPos){
    return findCount(altReadCount, targetPos);
}

std::vector<std::string> SubEdge::showEdge(std::string message){
    std::vector<std::string> result;
    for(CountTable::iterator edgeIter = refReadCount.begin() ; edgeIter != refReadCount.end() ; edgeIter++ ){
        result.push_back(message +" -> ref_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second) + "];");
    }
    for(CountTable::iterator edgeIter = altReadCount.begin() ; edgeIter != altReadCount.end() ; edgeIter++ ){
        result.push_back(message +" -> alt_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second) + "];");
    }
    return result;
}

std::vector<std::pair<int,int>> SubEdge::getConnectPos(){
    std::vector<std::pair<int,int>> result;
    for(CountTable::iterator edgeIter = refReadCount.begin() ; edgeIter != refReadCount.end() ; edgeIter++ ){
        result.push_back( std::make_pair( (*edgeIter).first, 0 ) );
    }
    for(CountTable::iterator edgeIter = altReadCount.begin() ; edgeIter != altReadCount.end() ; edgeIter++ ){
        result.push_back( std::make_pair( (*edgeIter).first, 1 ) );
    }
    return result;
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

// rr, ra, ar, aa: read weight from this variant's REF/ALT allele to the target's REF/ALT allele
std::pair<PosAllele,PosAllele> VairiantGraph::findBestEdgePair(int currPos, int targetPos, float rr, float ra, float ar, float aa, double edgeThreshold, bool debug, int currType, int targetType, VoteResult &vote){
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
    if((currType == 0 && targetType == 2)||(currType == 2 && targetType == 0)){
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
    if ( (edgeSimilarRatio <= 0.1 && (rr + aa + ra + ar) >= 1)  || ((rr+aa)<1&&(ra+ar)>=1) || ((rr+aa)>=1&&(ra+ar)<1) ) {
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
std::pair<float,float> VairiantGraph::Onelongcase( const std::vector<VoteResult> &vote ){

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
        else if ( vote[i].ESR < 0.2 && vote[i].weight >= 1 && variantType[variantIndex(vote[i].Pos)] != 3 ) {
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
void VairiantGraph::scanVariantsAndBuildBlocks(std::vector<int>& hpResult, PhasedBlocks& phasedBlocks, std::string& outDotResult){
    const int n = variantPositions.size();
    // votes received by each variant: < HP1, HP2 > weight and the votes themselves
    std::vector<std::array<float, 2> > hpCountMap2(n, std::array<float, 2>{0.0f, 0.0f});
    std::vector<std::vector<VoteResult> > hpCountMap3(n);
    hpResult.assign(n, 0);
    int blockStart = -1;
    int lastConnectPos = -1;
    // Unless --sv-join-blocks: a variant whose votes all come from SVs starts
    // a new block, so an SV is never the only link between two blocks. SVs
    // still vote inside a block and are phased in it.
    std::vector<char> nonSvVote(params->svJoinBlocks ? 0 : n, 0);
    int lastConnectPosNonSv = -1;

    for(int curr = 0; curr + 1 < n; curr++ ){
        int next = curr + 1;
        int currPos = variantPositions[curr];
        int nextPos = variantPositions[next];

        // votes for earlier positions are never read again
        if(curr > 0)
            std::vector<VoteResult>().swap(hpCountMap3[curr - 1]);

        if(std::abs(nextPos - currPos) > params->distance){
            continue;
        }

        float h1 = hpCountMap2[curr][0];
        float h2 = hpCountMap2[curr][1];

        std::pair<float, float> special = Onelongcase(hpCountMap3[curr]);
        if(special.first != -1){
            h1 = special.first;
            h2 = special.second;
        }

        // Record per-variant vote diagnostics (h1, h2, Shannon entropy).
        // Later consumed in exportResult() and written to VCF as INFO
        // fields H1 / H2 / PE. Only positions that eventually get phased
        // will actually reach the VCF; this map may hold values for
        // positions that don't, and that's harmless.
        {
            float total = h1 + h2;
            float entropy = 0.0f;
            if(total > 0.0f){
                float p1 = h1 / total;
                float p2 = h2 / total;
                if(p1 > 0.0f) entropy -= p1 * std::log2(p1);
                if(p2 > 0.0f) entropy -= p2 * std::log2(p2);
            }
            h1weight[curr] = h1;
            h2weight[curr] = h2;
            variantEntropy[curr] = entropy;
            hasEntropy[curr] = 1;
        }

        bool onlySvVotes = !params->svJoinBlocks && h1 != h2 && !nonSvVote[curr] &&
                           currPos >= lastConnectPosNonSv;
        if(h1 == h2 || onlySvVotes){
            if(currPos < lastConnectPos && !onlySvVotes){
                continue;
            }

            blockStart = currPos;
            phasedBlocks[blockStart].push_back(curr);
            hpResult[curr] = 1;
        }
        else{
            int currHP = (h1 > h2 ? 1 : 2);
            hpResult[curr] = currHP;
            phasedBlocks[blockStart].push_back(curr);
        }

        if(!hasEdge[curr]){
            continue;
        }

        // weights of curr's GNN edges, converted once per value
        float gnnRaw[2] = {0.0f, 0.0f}, gnnWeight[2] = {0.0f, 0.0f};
        bool gnnWeightSet = false;

        for(int i = 0; i < params->connectAdjacent; i++ ){
            VoteResult vote(currPos, 1);
            const std::array<float, 4> &counts = edgeCounts[(size_t)curr * params->connectAdjacent + i];
            std::pair<PosAllele,PosAllele> bestEdgePair = findBestEdgePair(currPos, variantPositions[next], counts[0], counts[1], counts[2], counts[3], params->edgeThreshold, false, variantType[curr], variantType[next], vote);

            if(variantType[curr] == 4){
                vote.weight = 0.1;
            }

            if(bestEdgePair.first.second != -1){
                recordVoteForNextPosition(
                    hpResult[curr],
                    next,
                    bestEdgePair.first.second,
                    vote,
                    hpCountMap2,
                    hpCountMap3);
                // variant type 1: SV
                if(!params->svJoinBlocks && variantType[curr] != 1){
                    nonSvVote[next] = 1;
                    lastConnectPosNonSv = variantPositions[next];
                }

                if(params->generateDot){
                    appendDotEdges(outDotResult, currPos, bestEdgePair,
                                   hpCountMap2[curr][0],
                                   hpCountMap2[curr][1]);
                }
                if(params->enableGNN){
                    // the same edges as the DOT file: curr.1 -> next.a and
                    // curr.2 -> next.(3-a)
                    const float h1 = hpCountMap2[curr][0], h2 = hpCountMap2[curr][1];
                    if(!gnnWeightSet || h1 != gnnRaw[0] || h2 != gnnRaw[1]){
                        gnnRaw[0] = h1;
                        gnnRaw[1] = h2;
                        gnnWeight[0] = gnnEdgeWeight(h1);
                        gnnWeight[1] = gnnEdgeWeight(h2);
                        gnnWeightSet = true;
                    }
                    gnnEdges.add(currPos, gnnWeight[0], gnnWeight[1],
                                 bestEdgePair.first.first, bestEdgePair.first.second);
                }

                lastConnectPos = variantPositions[next];
            }

            next++;
            if(next == n){
                break;
            }
        }
    }

    // The main loop exits before processing the last position (it needs
    // a next-position to iterate). If the last position received votes
    // from earlier neighbours, record its entropy too so it also gets
    // H1 / H2 / PE in the output VCF.
    if(n > 0){
        int last = n - 1;
        if(!hasEntropy[last]){
            float lh1 = hpCountMap2[last][0];
            float lh2 = hpCountMap2[last][1];
            float ltotal = lh1 + lh2;
            float lentropy = 0.0f;
            if(ltotal > 0.0f){
                float lp1 = lh1 / ltotal;
                float lp2 = lh2 / ltotal;
                if(lp1 > 0.0f) lentropy -= lp1 * std::log2(lp1);
                if(lp2 > 0.0f) lentropy -= lp2 * std::log2(lp2);
            }
            h1weight[last] = lh1;
            h2weight[last] = lh2;
            variantEntropy[last] = lentropy;
            hasEntropy[last] = 1;
        }
    }
}

void VairiantGraph::materializeBlockResults(
    const std::vector<int>& hpResult,
    const PhasedBlocks& phasedBlocks){
    for(auto blockIter = phasedBlocks.begin(); blockIter != phasedBlocks.end(); blockIter++ ){
        if((*blockIter).second.size() <= 1){
            continue;
        }

        for(auto currIter = (*blockIter).second.begin(); currIter != (*blockIter).second.end(); currIter++ ){
            auto nextIter = std::next(currIter, 1);
            if(nextIter == (*blockIter).second.end()){
                continue;
            }

            // allele slots: variant index * 2 + (0 ref, 1 alt)
            int refStart = (*currIter) * 2;
            int altStart = refStart + 1;
            int refEnd = (*nextIter) * 2;
            int altEnd = refEnd + 1;

            bkResult[refStart] = (*blockIter).first + 1;
            bkResult[refEnd]   = (*blockIter).first + 1;
            bkResult[altStart] = (*blockIter).first + 1;
            bkResult[altEnd]   = (*blockIter).first + 1;

            if(currIter == (*blockIter).second.begin()){
                subNodeHP[refStart] = 0;
                subNodeHP[altStart] = 1;
            }

            if(hpResult[*currIter] == 0 || hpResult[*nextIter] == 0){

            }
            else if(hpResult[*currIter] == hpResult[*nextIter]){
                subNodeHP[refEnd] = subNodeHP[refStart];
                subNodeHP[altEnd] = subNodeHP[altStart];
            }
            else{
                subNodeHP[refEnd] = subNodeHP[altStart];
                subNodeHP[altEnd] = subNodeHP[refStart];
            }
        }
    }
}

void VairiantGraph::edgeConnectResult(){
    std::vector<int> hpResult;
    PhasedBlocks phasedBlocks;

    scanVariantsAndBuildBlocks(hpResult, phasedBlocks, dotResult);
    materializeBlockResults(hpResult, phasedBlocks);
}

int VairiantGraph::variantIndex(int position) const{
    auto it = std::lower_bound(variantPositions.begin(), variantPositions.end(), position);
    if( it == variantPositions.end() || *it != position )
        return -1;
    return it - variantPositions.begin();
}

VairiantGraph::VairiantGraph(PhasingParameters &in_params, std::string &in_chrName){
    params=&in_params;
    chrName = &in_chrName;
}

VairiantGraph::~VairiantGraph(){
}

void VairiantGraph::destroy(){
    dotResult.clear();
    dotResult.shrink_to_fit();
    GnnEdgeList().swap(gnnEdges);

    std::vector<std::array<float, 4> >().swap(edgeCounts);
    std::vector<char>().swap(hasEdge);
    std::vector<int>().swap(variantPositions);
    std::vector<int8_t>().swap(variantType);
    std::vector<int>().swap(bkResult);
    std::vector<int8_t>().swap(subNodeHP);
    std::vector<float>().swap(h1weight);
    std::vector<float>().swap(h2weight);
    std::vector<float>().swap(variantEntropy);
    std::vector<char>().swap(hasEntropy);
}

//check if the position is in the range of the cnv
bool VairiantGraph::isPositionInRange(int position, int start, int end){
    return position >= start && position <= end;
}

//calculate the mismatch rate of the cnv
void VairiantGraph::calculateCnvMismatchRate(const std::vector<ReadVariant>& in_readVariant, const Clip &clip, ReadCnvMismatch& readMmrate){
    //if the read variant is empty or the cnv is empty, return
    if(in_readVariant.empty() || clip.cnvVec.empty()){
        return;
    }
    size_t cnvIndex = 0;

    for(size_t r = 0; r < in_readVariant.size(); r++){
        const ReadVariant& read = in_readVariant[r];
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
                //if the variant position is in the range of the cnv and the allele is reference, increment the mismatch rate.
                //SV alleles (quality -1) are not mismatches: the clips of an SV's own
                //split alignments make the CNV interval at the SV (issue #126)
                if(isPositionInRange(variant.position, clip.cnvVec[i].first, clip.cnvVec[i].second) && variant.allele == 1 &&
                   variant.quality != -1){
                    readMmrate[r][clip.cnvVec[i].first]++;
                }
            }
            i++;
        }
        //update the index of the cnv vector
        cnvIndex = i > 0 ? i - 1 : 0;
    }
    
}

//aggregate the mismatch rate of the cnv
void VairiantGraph::aggregateCnvReadMismatchRate(const std::vector<ReadVariant>& in_readVariant, const Clip &clip, const ReadCnvMismatch& readMmrate, std::map<int, std::map<int, std::vector<int>>>& cnvReadMmrate) { 
    //if the read variant is empty or the cnv is empty, return
    if(in_readVariant.empty() || clip.cnvVec.empty()){
        return;
    }

    size_t cnvIndex = 0;
    const std::map<int, int> noMmrate;

    for(size_t r = 0; r < in_readVariant.size(); r++){
        const ReadVariant& read = in_readVariant[r];
        //if the read has no variant, continue
        if(read.variantVec.empty()){
            continue;
        }
        auto mmIter = readMmrate.find(r);
        const std::map<int, int>& cnv_mmrate_map = mmIter != readMmrate.end() ? mmIter->second : noMmrate;

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
                if(isPositionInRange(variant.position, clip.cnvVec[i].first, clip.cnvVec[i].second) && cnv_mmrate_map.find(clip.cnvVec[i].first) != cnv_mmrate_map.end()){
                    cnvReadMmrate[variant.position][variant.allele].push_back(cnv_mmrate_map.at(clip.cnvVec[i].first));
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
                    //if the miss rate is greater than 0.7, erase the variant.
                    //SV alleles (quality -1) are kept: the clips of an SV's own
                    //split alignments end a CNV interval at the SV, and its
                    //carriers count the SV itself as a mismatch (issue #126)
                    if(missIter != missRateMap.end() && missIter->second >= 0.7 &&
                       variantIter->quality != -1){
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
    
// Indices of the reads grouped by read_id (so in read name order), each read's
// alignments in their input order.
static std::vector<uint32_t> alignmentsByRead(const std::vector<ReadVariant> &in_readVariant){
    uint32_t maxId = 0;
    for( const ReadVariant &read : in_readVariant ) maxId = std::max(maxId, read.read_id);
    std::vector<uint32_t> next(in_readVariant.empty() ? 0 : (size_t)maxId + 1, 0);
    for( const ReadVariant &read : in_readVariant ) next[read.read_id]++;
    uint32_t sum = 0;
    for( uint32_t &c : next ){ uint32_t k = c; c = sum; sum += k; }
    std::vector<uint32_t> order(in_readVariant.size());
    for( size_t i = 0; i < in_readVariant.size(); i++ ) order[next[in_readVariant[i].read_id]++] = i;
    return order;
}

void VairiantGraph::filterOverlappingAlignments(std::vector<ReadVariant> &in_readVariant){
    // record need del read index
    std::vector<int> delReadIdx;
    // the kept alignments of the current read
    std::vector<int> readIdxVecRef;

    // Check for overlaps among different alignments of a read and filter out the shorter overlapping alignments.
    // One read at a time, its alignments in input order.
    std::vector<uint32_t> order = alignmentsByRead(in_readVariant);
    for (size_t g = 0; g < order.size(); ) {
      const uint32_t readId = in_readVariant[order[g]].read_id;
      // the first and last variant position of the read; {0, 0} before its
      // first alignment, so the overlap condition below does not hold then
      std::pair<int,int> readRange(0, 0);
      readIdxVecRef.clear();
      for ( ; g < order.size() && in_readVariant[order[g]].read_id == readId; g++) {
        int readIter = order[g];
        int is_toDelete = 0;
        int firstVariantPos = in_readVariant[readIter].variantVec.front().position;
        int lastVariantPos = in_readVariant[readIter].variantVec.back().position;

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
    }

    // sort read index
    std::sort(delReadIdx.begin(), delReadIdx.end());
    // remove overlap alignment
    delReadIdx.push_back((int)in_readVariant.size());
    int saveIter = *(delReadIdx.begin());
    for (auto delIter = delReadIdx.begin(), nextdelIter = std::next(delReadIdx.begin(), 1); nextdelIter != delReadIdx.end(); delIter++ , nextdelIter++) {
        auto nowDelIter = *delIter+1;
        while (nowDelIter<*nextdelIter){
            in_readVariant[saveIter++]=std::move(in_readVariant[nowDelIter++]);
        }
    }
    in_readVariant.erase( std::next(in_readVariant.begin(), saveIter), in_readVariant.end()); 
}

void VairiantGraph::applyCnvFilter(std::vector<ReadVariant> &in_readVariant, Clip &clip){
    CnvStatistics cnvStats;
    //calculate the mismatch rate of the cnv
    ReadCnvMismatch readMmrate;
    calculateCnvMismatchRate(in_readVariant, clip, readMmrate);
    //aggregate the mismatch rate of the cnv
    aggregateCnvReadMismatchRate(in_readVariant, clip, readMmrate, cnvStats.cnvReadMmrate);
    ReadCnvMismatch().swap(readMmrate);
    //calculate the average mismatch rate of the cnv
    calculateAverageMismatchRate(clip, cnvStats.cnvReadMmrate, cnvStats.missRateMap);
    //filter the variants with high mismatch rate
    filterHighMismatchVariants(in_readVariant, clip, cnvStats.missRateMap);
}

// SV, modification and indel calls carry their type in the quality field; return the type
// (0 SNP, 1 SV, 2 modification, 3 indel, 4 danger indel) and set the edge quality
static int normalizeVariantQuality(Variant &variant){
    if( variant.quality == -2 || variant.quality == -3 ){ variant.quality = 60; return 2; }
    if( variant.quality == -1 ){ variant.quality = (variant.allele == 1 ? 60 : 30); return 1; }
    if( variant.quality == -4 ){ variant.quality = 60; return 3; }
    if( variant.quality == -5 ){ variant.quality = 60; return 4; }
    return 0;
}

void VairiantGraph::buildVariantGraph(std::vector<ReadVariant> &in_readVariant){
    // variant positions in genomic order
    variantPositions.clear();
    for( const ReadVariant &read : in_readVariant ){
        for( const Variant &variant : read.variantVec ){
            variantPositions.push_back(variant.position);
        }
    }
    std::sort(variantPositions.begin(), variantPositions.end());
    variantPositions.erase(std::unique(variantPositions.begin(), variantPositions.end()), variantPositions.end());
    variantPositions.shrink_to_fit();

    const size_t n = variantPositions.size();
    edgeCounts.assign(n * params->connectAdjacent, std::array<float, 4>{0.0f, 0.0f, 0.0f, 0.0f});
    hasEdge.assign(n, 0);
    variantType.assign(n, 0);
    bkResult.assign(n * 2, 0);
    subNodeHP.assign(n * 2, 0);
    h1weight.assign(n, 0.0f);
    h2weight.assign(n, 0.0f);
    variantEntropy.assign(n, 0.0f);
    hasEntropy.assign(n, 0);

    // variant types (the last read seen sets the type of a position)
    for( const ReadVariant &read : in_readVariant ){
        for( auto variant : read.variantVec ){
            variantType[variantIndex(variant.position)] = normalizeVariantQuality(variant);
        }
    }

    // alignments of one read, in read name order
    std::vector<uint32_t> order = alignmentsByRead(in_readVariant);

    ReadVariant merged;
    for(size_t g = 0; g < order.size(); ){
        const uint32_t readId = in_readVariant[order[g]].read_id;
        merged.variantVec.clear();
        size_t e = g;
        for( ; e < order.size() && in_readVariant[order[e]].read_id == readId; e++){
            for( auto variant : in_readVariant[order[e]].variantVec ){
                normalizeVariantQuality(variant);
                merged.variantVec.push_back(variant);
            }
        }
        g = e;
        // the map only held reads with variants
        if( merged.variantVec.empty() )
            continue;
        ReadVariant *readIter = &merged;
        readIter->sort();
        
        // iter all pair of snp and construct initial graph
        std::vector<Variant>::iterator variant1Iter = readIter->variantVec.begin();
        std::vector<Variant>::iterator variant2Iter = std::next(variant1Iter,1);
        // index of variant1Iter; the read's variants are sorted, so the search starts there
        int index1 = variantIndex(variant1Iter->position);

        while(variant1Iter != readIter->variantVec.end() && variant2Iter != readIter->variantVec.end() ){
            while( variantPositions[index1] < variant1Iter->position ) index1++;
            hasEdge[index1] = 1;
            int index2 = index1;

            // add edge process
            for(int nextNode = 0 ; nextNode < params->connectAdjacent; nextNode++){
                while( variantPositions[index2] < variant2Iter->position ) index2++;
                int distance = index2 - index1;
                // only the next connectAdjacent variants are ever compared with this one
                if( distance >= 1 && distance <= params->connectAdjacent
                    && (variant1Iter->allele == 0 || variant1Iter->allele == 1)
                    && (variant2Iter->allele == 0 || variant2Iter->allele == 1) ){
                    float &count = edgeCounts[(size_t)index1 * params->connectAdjacent + distance - 1][variant1Iter->allele * 2 + variant2Iter->allele];
                    //if the base quality on both snps is high enough, the edge has normal weight
                    if ( variant1Iter->quality >= params->baseQuality && variant2Iter->quality >= params->baseQuality)
                        count++;
                    else
                        count = count + params->edgeWeight ;
                }

                // next snp
                variant2Iter++;
                if( variant2Iter == readIter->variantVec.end() ){
                    break;
                }
            }

            variant1Iter++;
            variant2Iter = std::next(variant1Iter,1);
        }

        //count the ref and alt base amount of the last variant on the read
        if ( variant1Iter != readIter->variantVec.end() && variant2Iter == readIter->variantVec.end() ) {
            while( variantPositions[index1] < variant1Iter->position ) index1++;
            hasEdge[index1] = 1;
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
    int index = variantIndex(variant.position);
    if(index < 0 || (variant.allele != 0 && variant.allele != 1)){
        return;
    }
    // allele slot of this read's allele
    int refAllele = index * 2 + variant.allele;

    if(bkResult[refAllele] == 0){
        return;
    }

    int refHaplotype = subNodeHP[refAllele];
    int positionVariantType = variantType[index];

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

            for(auto variantIter = (*readIter).variantVec.begin(); variantIter != (*readIter).variantVec.end(); variantIter++){
                if((*variantIter).allele == 0 || (*variantIter).allele == 1){
                    int index = variantIndex((*variantIter).position);
                    if(index >= 0)
                        hpAlleleCountMap[belongHP][index][(*variantIter).allele]++;
                }
            }
        }
    }
}

void VairiantGraph::reassignVariantHaplotypes(const HpAlleleCountMap& hpAlleleCountMap){
    auto getHpAlleleCount = [&](int haplotype, int index, int allele){
        return hpAlleleCountMap[haplotype][index][allele];
    };

    double snpConfidenceThreshold = params->snpConfidence;

    std::fill(subNodeHP.begin(), subNodeHP.end(), 0);

    // reassign allele result
    for(int index = 0; index < (int)variantPositions.size(); index++){
        int refAllele = index * 2;
        int altAllele = refAllele + 1;

        double hp1Ref = getHpAlleleCount(0, index, 0);
        double hp1Alt = getHpAlleleCount(0, index, 1);
        double hp2Ref = getHpAlleleCount(1, index, 0);
        double hp2Alt = getHpAlleleCount(1, index, 1);
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
            subNodeHP[refAllele] = hp1Result;
            subNodeHP[altAllele] = hp2Result;
        }
        else{
            bkResult[refAllele] = 0;
            bkResult[altAllele] = 0;
        }
    }
}

void VairiantGraph::readCorrection(){
    HpAlleleCountMap hpAlleleCountMap;
    for(auto &counts : hpAlleleCountMap)
        counts.assign(variantPositions.size(), std::array<double, 2>{0.0, 0.0});
    assignReadHaplotypesAndCollectAlleleCounts(hpAlleleCountMap);
    reassignVariantHaplotypes(hpAlleleCountMap);
}

void VairiantGraph::writingDotFile(std::string dotPrefix){

    // Prefix dot files with the run's --out-prefix, e.g.
    //   result.chr1.dot   instead of just   chr1.dot
    std::string outPath = params->resultPrefix + "." + dotPrefix + ".dot";
    std::ofstream resultVcf(outPath);

    if(!resultVcf.is_open()){
        std::cerr<< "Fail to open write file: " << outPath << "\n";
    }
    else{
        resultVcf << "digraph G {\n";

        resultVcf << dotResult;
        resultVcf << "}\n";
    }
    return;
}

void VairiantGraph::exportGnnEdges(GnnEdgeList &out){
    // drop the spare capacity; these are kept until the GNN runs
    gnnEdges.shrink_to_fit();
    out.swap(gnnEdges);
}

void VairiantGraph::exportResult(std::string chrName, PhasingResult &result){
    
    // loop all position
    for( int index = 0 ; index < (int)variantPositions.size() ; index++ ){
        int refAllele = index * 2;
        int altAllele = refAllele + 1;

        // 0: not in a block
        int block = (bkResult[refAllele] != 0 ? bkResult[refAllele] : bkResult[altAllele]);
        if( block == 0 )
            continue;

        PhasingElement tmp;
        tmp.block = block;
        tmp.RAstatus = std::to_string(subNodeHP[refAllele]) + "|" + std::to_string(subNodeHP[altAllele]);

        // Attach vote diagnostics (h1, h2, entropy) so writeDataLine
        // can emit INFO/H1, INFO/H2, INFO/PE for phased positions.
        // These are filled by scanVariantsAndBuildBlocks().
        tmp.h1      = h1weight[index];
        tmp.h2      = h2weight[index];
        tmp.entropy = variantEntropy[index];

        std::string key = chrName + "_" + std::to_string( variantPositions[index] );
        result[key] = tmp;
    }
}

int VairiantGraph::totalNode(){
    return variantPositions.size();
}

void VairiantGraph::phasingProcess(){
    // This step involves converting all reads into a graph structure, which will be stored as an edge list
    // in a two-layer map. The first layer of the map uses the starting coordinate as the key and contains
    // a second layer map as the value. The second layer map uses the destination coordinate as the key and
    // stores the number of support read as values. (There is another map used for debugging purposes that
    // treats the read name vector as a value.) The method begins by visiting the coordinates covered by each 
    // read and recording this information in 'variantPositions'. Subsequently, it connects the coordinates contained 
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
