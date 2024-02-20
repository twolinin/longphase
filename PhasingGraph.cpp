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
	if ( currentQuality >= baseQuality && connectNode.quality >= baseQuality )
            (*refReadCount)[connectNode.position]++;
        else {
            (*refReadCount)[connectNode.position] = (*refReadCount)[connectNode.position] + edgeWeight ;
        }
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
	if ( currentQuality >= baseQuality && connectNode.quality >= baseQuality )
            (*altReadCount)[connectNode.position]++;
        else {
            (*altReadCount)[connectNode.position] = (*altReadCount)[connectNode.position] + edgeWeight ;
        }
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

VariantEdge::VariantEdge(int inCurrPos){
    currPos = inCurrPos;
    alt = new SubEdge();
    ref = new SubEdge();
}

//VariantEdge
std::pair<PosAllele,PosAllele> VariantEdge::findBestEdgePair(int targetPos, bool isONT, double edgeThreshold, bool debug){
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
    
    if( edgeSimilarRatio > edgeThreshold ){
        refAllele = -1;
        altAllele = -1;
    }
    
    if(debug){
        std::cout<< currPos << "\t->\t" << targetPos << "\t|rr aa | ra ar\t" << "\t" << rr << "\t" << aa << "\t" << ra << "\t" << ar  << "\n";
    }
    
    // create edge pairs
    PosAllele refEdge = std::make_pair( targetPos, refAllele );
    PosAllele altEdge = std::make_pair( targetPos, altAllele );
    // return edge pair
    return std::make_pair( refEdge, altEdge );
}

std::pair<int,int> VariantEdge::findNumberOfRead(int targetPos){
    std::pair<int,int> refBestPair  = ref->BestPair(targetPos);
    std::pair<int,int> altBestPair  = alt->BestPair(targetPos);
    // get the weight of each pair
    int rr = refBestPair.first;
    int ra = refBestPair.second;
    int ar = altBestPair.first;
    int aa = altBestPair.second;
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

//VairiantGraph
void VairiantGraph::edgeConnectResult(){
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
        int h1 = (*hpCountMap)[currPos][1].size();
        int h2 = (*hpCountMap)[currPos][2].size();

        // new block, set this position as block start 
        if( h1 == 0 && h2 == 0 ){
            // No new blocks should be created if the next SNP has already been picked up
            if( currPos < lastConnectPos ){
                continue;
            }
            
            blockStart = currPos;
            (*phasedBlocks)[blockStart].push_back(currPos);
            (*hpResult)[currPos] = 1;
        }
        else{
            if( h1 > h2 || h1 < h2 ){
                int currHP = ( h1 > h2 ? 1 : 2 );
                (*hpResult)[currPos] = currHP;
                (*phasedBlocks)[blockStart].push_back(currPos);
            }
        }
        
        // Check if there is no edge from current node
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList->find( currPos );
        if( edgeIter==edgeList->end() ){
            continue;
        }
        
        // check connect between surrent SNP and next n SNPs
        for(int i = 0 ; i < params->connectAdjacent ; i++ ){
            // consider reads from the currnt SNP and the next (i+1)'s SNP
            std::pair<PosAllele,PosAllele> tmp = edgeIter->second->findBestEdgePair(nextNodeIter->first, params->isONT, params->edgeThreshold, false);
            // -1 : no connect  
            //  1 : the haplotype of next (i+1)'s SNP are same as previous
            //  2 : the haplotype of next (i+1)'s SNP are different as previous
            if( tmp.first.second != -1 ){
                // record the haplotype resut of next (i+1)'s SNP
                if( (*hpResult)[currPos] == 1 ){
                    if( tmp.first.second == 1 ){
                        (*hpCountMap)[nextNodeIter->first][1].push_back(currPos);
                    }
                    if( tmp.first.second == 2 ){
                        (*hpCountMap)[nextNodeIter->first][2].push_back(currPos);
                    }
                }
                if( (*hpResult)[currPos]==2 ){
                    if( tmp.first.second == 1 ){
                        (*hpCountMap)[nextNodeIter->first][2].push_back(currPos);
                    }
                    if( tmp.first.second == 2 ){
                        (*hpCountMap)[nextNodeIter->first][1].push_back(currPos);
                    }
                }
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

void VairiantGraph::addEdge(std::vector<ReadVariant> &in_readVariant){

    readVariant = &in_readVariant;
    // iter all read
    for(std::vector<ReadVariant>::iterator readIter = in_readVariant.begin() ; readIter != in_readVariant.end() ; readIter++ ){
        // Creating a pseudo read which allows filtering out variants that should not be phased
        ReadVariant tmpRead;
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
            // The remaining variants will be labeled as SNPs
            else{
                (*variantType)[variant.position] = 0;
            }
            
            tmpRead.variantVec.push_back(variant);
            
            // Each position will record the included reads and their corresponding base qualities.
            auto variantIter = totalVariantInfo->find(variant.position);
            
            if( variantIter == totalVariantInfo->end() ){
                (*totalVariantInfo)[variant.position] = new ReadBaseMap();
            }
            
            (*(*totalVariantInfo)[variant.position])[(*readIter).read_name] = variant.quality;
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
                    (*edgeList)[(*variant1Iter).position]->ref->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).read_name,params->baseQuality,params->edgeWeight);
                // this allele support alt
                if( (*variant1Iter).allele == 1 )
                    (*edgeList)[(*variant1Iter).position]->alt->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).read_name,params->baseQuality,params->edgeWeight);
                
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

void VairiantGraph::readCorrection(){
    // haplotype, <position <allele, base count>>
    std::map<int,std::map<int,std::map<int,int>>> *hpAlleleCountMap = new std::map<int,std::map<int,std::map<int,int>>>;

    // iter all read, determine the haplotype of the read
    for(std::vector<ReadVariant>::iterator readIter = (*readVariant).begin() ; readIter != (*readVariant).end() ; readIter++ ){
        double refCount = 0;
        double altCount = 0;
        
        // loop all variant 
        for( auto variant : (*readIter).variantVec ){
            PosAllele refAllele = std::make_pair( variant.position , variant.allele+1);
            std::map<PosAllele,int>::iterator nodePS = bkResult->find(refAllele);

            if( nodePS != bkResult->end() ){
                if( (*bkResult)[refAllele] != 0 ){
                    if((*subNodeHP)[refAllele]==0)refCount++;
                    else altCount++;
                }
            }
        }
        
        // tag high confident reads
        if( std::max(refCount,altCount)/(refCount+altCount) > params->readConfidence ){
            // tag read with the corresponding haplotype
            int belongHP = ( refCount > altCount ? 0 : 1 );
            (*readHpMap)[(*readIter).read_name] = belongHP;
            
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


