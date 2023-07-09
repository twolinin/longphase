#include "PhasingGraph.h"

//SubEdge

SubEdge::SubEdge():readCount(0){ 
    refRead = new std::map<int, std::vector<std::string> >;
    altRead = new std::map<int, std::vector<std::string> >;
    refQuality = new std::map<int, int>;
    altQuality = new std::map<int, int>;
    refReadCount = new std::map<int, int>;
    altReadCount = new std::map<int, int>;
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

void SubEdge::addSubEdge(int currentQuality, Variant connectNode, std::string readName){
    if(connectNode.allele == 0 ){
        (*refRead)[connectNode.position].push_back(readName);
        // quality sum
        std::map<int, int>::iterator rqIter = refQuality->find(connectNode.position);
        if( rqIter == refQuality->end() ){
            (*refQuality)[connectNode.position] = currentQuality + connectNode.quality;
            (*refReadCount)[connectNode.position] = 1;
        }
        else{
            (*refQuality)[connectNode.position] += currentQuality + connectNode.quality;
            (*refReadCount)[connectNode.position]++;
        }
    }
    if(connectNode.allele == 1 ){
        (*altRead)[connectNode.position].push_back(readName);
        // quality sum
        std::map<int, int>::iterator aqIter = altQuality->find(connectNode.position);
        if( aqIter == altQuality->end() ){
            (*altQuality)[connectNode.position] = currentQuality + connectNode.quality;
            (*altReadCount)[connectNode.position] = 1;
        }
        else{
            (*altQuality)[connectNode.position] += currentQuality + connectNode.quality;
            (*altReadCount)[connectNode.position]++;
        }
    }
    readCount++;
}

std::vector<std::string> SubEdge::showEdge(std::string message){
    std::vector<std::string> result;
    for(std::map<int, std::vector<std::string> >::iterator edgeIter = refRead->begin() ; edgeIter != refRead->end() ; edgeIter++ ){
        result.push_back(message +" -> ref_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second.size()) + "];");
    }
    for(std::map<int, std::vector<std::string> >::iterator edgeIter = altRead->begin() ; edgeIter != altRead->end() ; edgeIter++ ){
        result.push_back(message +" -> alt_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second.size()) + "];");
    }
    return result;
}

std::vector<std::pair<int,int>> SubEdge::getConnectPos(){
    std::vector<std::pair<int,int>> result;
    for(std::map<int, std::vector<std::string> >::iterator edgeIter = refRead->begin() ; edgeIter != refRead->end() ; edgeIter++ ){
        result.push_back( std::make_pair( (*edgeIter).first, 0 ) );
    }
    for(std::map<int, std::vector<std::string> >::iterator edgeIter = altRead->begin() ; edgeIter != altRead->end() ; edgeIter++ ){
        result.push_back( std::make_pair( (*edgeIter).first, 1 ) );
    }
    return result;
}

std::pair<int,int> SubEdge::BestPair(int targetPos){
    // default, no read support target position
    int refResult = 0;
    int altResult = 0;
    // find the edge between local and target position 
    std::map<int, std::vector<std::string> >::iterator refEdgeIter = refRead->find(targetPos);
    std::map<int, std::vector<std::string> >::iterator altEdgeIter = altRead->find(targetPos);
    // local connect target REF allele
    if( refEdgeIter != refRead->end() ){
        refResult = (*refEdgeIter).second.size();
    }
    // local connect target ALT allele
    if( altEdgeIter != altRead->end() ){
        altResult = (*altEdgeIter).second.size();   
    }            
    return std::make_pair( refResult, altResult );
}

void SubEdge::getReadVariant(std::map< std::string,std::map<int,int> > &readVariantMap){

    for(std::map<int, std::vector<std::string> >::iterator edgeIter = refRead->begin() ; edgeIter != refRead->end() ; edgeIter++ ){
        for(std::vector<std::string>::iterator readIter = (*edgeIter).second.begin() ; readIter != (*edgeIter).second.end() ; readIter++ ){
            readVariantMap[(*readIter)][(*edgeIter).first]=0;
        }
    }
    
    for(std::map<int, std::vector<std::string> >::iterator edgeIter = altRead->begin() ; edgeIter != altRead->end() ; edgeIter++ ){
        for(std::vector<std::string>::iterator readIter = (*edgeIter).second.begin() ; readIter != (*edgeIter).second.end() ; readIter++ ){
            readVariantMap[(*readIter)][(*edgeIter).first]=1;
        }
    }
}

void SubEdge::deletEdge(int pos){
	for(std::map<int, std::vector<std::string> >::iterator edgeIter = refRead->begin() ; edgeIter != refRead->end() ; edgeIter++ ){
        //result.push_back(message +" -> ref_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second.size()) + "];");
		if( (*edgeIter).first - pos > 40000 && (*edgeIter).second.size()<=1){
			std::cout<<"deletEdge\t"<<pos<<"\t"<<(*edgeIter).first<<"\t"<<(*edgeIter).second.size()<<"\n";
			(*edgeIter).second.clear();
			//(*refReadIter).second.shrink_to_fit();
		}
    }
    for(std::map<int, std::vector<std::string> >::iterator edgeIter = altRead->begin() ; edgeIter != altRead->end() ; edgeIter++ ){
        //result.push_back(message +" -> alt_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second.size()) + "];");
		if( (*edgeIter).first - pos > 40000 && (*edgeIter).second.size()<=1){
			std::cout<<"deletEdge\t"<<pos<<"\t"<<(*edgeIter).first<<"\t"<<(*edgeIter).second.size()<<"\n";
			(*edgeIter).second.clear();
			//(*refReadIter).second.shrink_to_fit();
		}
	}
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

int SubEdge::getRefReadCount(int targetPos){
    std::map<int, int>::iterator posIter = refReadCount->find(targetPos);
    if( posIter != refReadCount->end() ){
        return (*refReadCount)[targetPos];
    }
    return 0;
}

int SubEdge::getAltReadCount(int targetPos){
    std::map<int, int>::iterator posIter = altReadCount->find(targetPos);
    if( posIter != altReadCount->end() ){
        return (*altReadCount)[targetPos];
    }
    return 0;
}

VariantEdge::VariantEdge(int inCurrPos){
    currPos = inCurrPos;
    alt = new SubEdge();
    ref = new SubEdge();
}

//VariantEdge
std::pair<PosAllele,PosAllele> VariantEdge::findBestEdgePair(int targetPos, bool isONT, double diffRatioThreshold, bool repair, bool debug){
    std::pair<int,int> refBestPair  = ref->BestPair(targetPos);
    std::pair<int,int> altBestPair  = alt->BestPair(targetPos);
    // get the weight of each pair
    int rr = refBestPair.first;
    int ra = refBestPair.second;
    int ar = altBestPair.first;
    int aa = altBestPair.second;
    // initialize the edge connection
    // -1 : not connect
    int refAllele = -1;
    int altAllele = -1;

    // If two choice have the same weight, 
    // the two edges have the weight is more reliable than one.
    // i.e. methylation
    if(isONT){
        if( rr+aa == ra +ar ){
            if(rr!=0 && aa !=0){
                rr+=1;
                aa+=1;
            }
            if(ra!=0 && ar !=0){
                ra+=1;
                ar+=1;
            }
        }
    }
    
    double readsThreshold = (double)std::min((rr+aa),(ar+ra)) / (double)std::max((rr+aa),(ar+ra));

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
    
    if( readsThreshold > diffRatioThreshold && repair){
        refAllele = -1;
        altAllele = -1;
    }
    
    if(debug){
        std::cout<< "rr aa | ra ar\t" << "\t" << rr << "\t" << aa << "\t" << ra << "\t" << ar  << "\n";
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
    // current snp, result haplotype (1 or 2)
    std::map<int,int> *inconsistentSnpMap = new std::map<int,int>;
    // < block start, <snp in this block> >
    std::map<int,std::vector<int> > *phasedBlocks = new std::map<int,std::vector<int> >;
    
    int blockStart = -1;
    int prevSNP = -1;
    int currPos = -1;
    int nextPos = -1;
    int lastConnectPos = -1;
    
    // Visit all position and assign SNPs to haplotype.
    // Avoid recording duplicate information,
    // only one of the two alleles needs to be used for each SNP
    for(std::map<int,ReadBaseMap*>::iterator nodeIter = nodeInfo->begin() ; nodeIter != nodeInfo->end() ; nodeIter++ ){
        // check next position
        std::map<int,ReadBaseMap*>::iterator nextNodeIter = std::next(nodeIter, 1);
        if( nextNodeIter == nodeInfo->end() ){
             break;
        }
        
        prevSNP = currPos;
        currPos = nodeIter->first;
        nextPos = nextNodeIter->first;
        
        // check distance between two position 
        if(std::abs(nextPos-currPos) > params->distance ){
            continue;
        }
        
        // get the number of HP1 and HP2 supported reference allele
        int h1 = (*hpCountMap)[currPos][1].size();
        int h2 = (*hpCountMap)[currPos][2].size();
        
        double phasedConfidence   = (double)std::max(h1,h2)/(double)(h1+h2);
        double haplotypeInconsistentRatio = (double)std::min(h1,h2)/(double)(h1+h2);
        
        // check and record inconsistent snp
        // using haplotype ratio to determine whether inconsistency occurs
        // i.g. HP1 : HP2 = 10 : 2
        // the SNPs support HP2 are inconsistent
        if( h1!=0 && h2!=0 && haplotypeInconsistentRatio <= params->judgeInconsistent ){
            // calculate h1 ratio
            double h1ratio = (double)h1/(double)(h1+h2);

            // support h1 reads are inconsistent
            if( h1ratio == haplotypeInconsistentRatio ){
                for( unsigned int i = 0 ; i < (*hpCountMap)[currPos][1].size() ; i++ ){
                    (*inconsistentSnpMap)[(*hpCountMap)[currPos][1][i]]++;
                }
            }
            // support h2 reads are inconsistent
            else{
                for( unsigned int i = 0 ; i < (*hpCountMap)[currPos][2].size() ; i++ ){
                    (*inconsistentSnpMap)[(*hpCountMap)[currPos][2][i]]++;
                }
            }
        }
        
        // new block, set this position as block start 
        if( h1 == 0 && h2 == 0 ){
            
            // No new blocks should be created if the next SNP has already been picked up
            if( currPos < lastConnectPos ){
                continue;
            }
            
            blockStart=currPos;
            (*phasedBlocks)[blockStart].push_back(currPos);
            (*hpResult)[currPos] = 1;
        }
        // poorly classified SNPs
        // votes for this SNP from the previous N SNPs are tied or similar
        else if( 1 - params->confidentHaplotype <= phasedConfidence  && phasedConfidence <= params->confidentHaplotype ){
            // check previos SNP in edge list
            std::map<int,VariantEdge*>::iterator edgeIter = edgeList->find( prevSNP );
            if( edgeIter != edgeList->end() ){
                // record none haplotype
                std::pair<PosAllele,PosAllele> preResult = (*edgeIter).second->findBestEdgePair(currPos, params->isONT, params->readsThreshold, true, false);  
                if( preResult.first.second == 1 || preResult.first.second == 2 ){
                    // Although the h1 and h2 votes of the previous N points are tied, 
                    // but this position is not tied with the previous SNP, 
                    // push this SNP but not assign haplotype, this SNP may classified by readCorrection()
                    //(*hpResult)[currPos] = 0;	
                    //(*phasedBlocks)[blockStart].push_back(currPos);	
                    //change:ConnetTieByPrev	
                    if( (*hpResult)[(*phasedBlocks)[blockStart].back()] == 1 ){	
                        if( preResult.first.second == 1 ){	
                            (*hpResult)[currPos] = 1 ;	
                        }	
                        if( preResult.first.second == 2 ){	
                            (*hpResult)[currPos] = 2 ;	
                        }	
                    }	
                    else if( (*hpResult)[(*phasedBlocks)[blockStart].back()] == 2 ) {	
                        if( preResult.first.second == 1 ){	
                            (*hpResult)[currPos] = 2 ;	
                        }	
                        if( preResult.first.second == 2 ){	
                            (*hpResult)[currPos] = 1 ;	
                        }	
                    }
                    (*phasedBlocks)[blockStart].push_back(currPos);
                }
            }
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
        if( edgeIter==edgeList->end() )
            continue;
        
        
        // check connect between surrent SNP and next n SNPs
        for(int i = 0 ; i < params->connectAdjacent ; i++ ){
            // consider reads from the currnt SNP and the next (i+1)'s SNP
            std::pair<PosAllele,PosAllele> tmp = edgeIter->second->findBestEdgePair(nextNodeIter->first, params->isONT, 1, false, false);
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
            if( nextNodeIter == nodeInfo->end() )
                break;
        }
    }

    // Correction inconsistent SNPs
    std::vector<std::string>::iterator edgeIter = dotResult.begin();
    for(std::map<int,int>::iterator snpIter = inconsistentSnpMap->begin();snpIter != inconsistentSnpMap->end();snpIter++){
        if( (*snpIter).second >= params->inconsistentThreshold ){
            // change haplotype result
            (*hpResult)[(*snpIter).first] = ((*hpResult)[(*snpIter).first] == 1 ? 2 : 1);
            // modify dot file
            if( params->generateDot ){
                // speedup for read dotResult vector
                std::vector<std::string>::iterator edgeTmpIter = edgeIter;
                bool start_replace = false;
                int count = 100;
                
                while( count >= 0 && edgeTmpIter != dotResult.end() && edgeIter != dotResult.end() ){
                    std::string originPattern1 = std::to_string((*snpIter).first+1) + ".1\t->\t";
                    std::string originPattern2 = std::to_string((*snpIter).first+1) + ".2\t->\t";
                    if( (*edgeTmpIter).find(originPattern1) != std::string::npos || 
                        (*edgeTmpIter).find(originPattern2) != std::string::npos){
                            
                        if( (*edgeTmpIter)[(*edgeTmpIter).length()-1] == '1' ){
                            (*edgeTmpIter)[(*edgeTmpIter).length()-1] = '2';
                        }
                        else{
                           (*edgeTmpIter)[(*edgeTmpIter).length()-1] = '1'; 
                        }    
                        
                        if(!start_replace){
                            start_replace=true;
                            edgeIter = edgeTmpIter;
                        }
                    }
                    
                    if(start_replace){
                        count--;
                    }
                    edgeTmpIter++;
                }
            }
        }
    }
    
    // loop all block
    for(auto blockIter = phasedBlocks->begin() ; blockIter != phasedBlocks->end() ; blockIter++ ){
        // check block size
        if( (*blockIter).second.size()<=1 ){
            continue;
        }
        
        // loop block's node
        for(auto currIter = (*blockIter).second.begin() ; currIter != (*blockIter).second.end() ; currIter++ ){
            // check next node
            auto nextIter = std::next(currIter,1);
            if( nextIter == (*blockIter).second.end() )
                continue;

            PosAllele refStart = std::make_pair((*currIter), 1);
            PosAllele refEnd;
            PosAllele altStart = std::make_pair((*currIter), 2);
            PosAllele altEnd;
            
            if( (*hpResult)[(*currIter)] == 0 || (*hpResult)[(*nextIter)] == 0 ){
                
            }
            else if( (*hpResult)[(*currIter)] == (*hpResult)[(*nextIter)] ){
                refEnd = std::make_pair((*nextIter), 1);
                altEnd = std::make_pair((*nextIter), 2);
            }
            else{
                refEnd = std::make_pair((*nextIter), 2);
                altEnd = std::make_pair((*nextIter), 1);
            }
            (*bestEdgeConnect)[refEnd] = refStart;
            (*bestEdgeConnect)[altEnd] = altStart;
        }
    }
    
    delete hpCountMap;
    delete hpResult;
    delete inconsistentSnpMap;
    delete phasedBlocks;
}

std::vector<ReadVariant> VairiantGraph::getBlockRead(std::pair<int,std::vector<int> > currentBlockVec, BlockRead &totalRead , int sampleNum){

    std::vector<ReadVariant> resultVector;
    // get ref reads
    std::map< std::string,std::map<int,int> > readVariantMap;
    
    // tag read belongs to which haplotype
    int SNPcount = 0;
    for(std::vector<int>::iterator posIter = currentBlockVec.second.begin(); posIter != currentBlockVec.second.end() ; posIter++ ){
        SNPcount++;
        // Because the read length cannot span all SNPs
        // a few SNPs before sampling are enough
        if( SNPcount > sampleNum )
            break;
            
        // prevent empty edge        
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList->find( (*posIter) );
        if( edgeIter==edgeList->end() )
            continue;
        
        (*edgeIter).second->ref->getReadVariant(readVariantMap);
        (*edgeIter).second->alt->getReadVariant(readVariantMap);

    }  
    
    //translate readVariantMap struct to ReadVariant struct
    for(auto readIter = readVariantMap.begin() ; readIter != readVariantMap.end() ; readIter++ ){
        ReadVariant tmpReadVariant;
        tmpReadVariant.read_name = (*readIter).first;
        for(auto variantIter = (*readIter).second.begin(); variantIter != (*readIter).second.end() ; variantIter++ ){
            Variant tmpVariant((*variantIter).first, (*variantIter).second, 30);
            tmpReadVariant.variantVec.push_back(tmpVariant);
        }
        resultVector.push_back(tmpReadVariant);
        
        totalRead.recordRead((*readIter).first);
    }
    return resultVector;
}
/*
bool VairiantGraph::blockPhaseing(){
    
    bool connect = false;
    
    // loop all block, prepare to create new edges between two block
    // PS number, < position vector >
    for(std::map<int,std::vector<int> >::iterator blockVecIter = blockVec->begin() ; blockVecIter != blockVec->end() ; blockVecIter++ ){
        // Number of SNPs checked
        //int sampleNum = 15;
        // get next block
        std::map<int,std::vector<int> >::iterator nextBlockIter = std::next(blockVecIter, 1);
        if( nextBlockIter == blockVec->end() )
            break;
        if( nextBlockIter == blockVecIter )
            break;
        // sort previous block
        std::sort((*blockVecIter).second.rbegin(), (*blockVecIter).second.rend());
        // sort following block
        std::sort((*nextBlockIter).second.begin(), (*nextBlockIter).second.end());
        // The last SNP of the previous block
        std::map<int,ReadBaseMap*>::iterator breakNode = nodeInfo->find((*blockVecIter).second[0]);
        // The last SNP of the following block
        std::map<int,ReadBaseMap*>::iterator nodeIter = nodeInfo->find((*nextBlockIter).second[0]);
        // The distance between the two blocks is abnormal, which may be affected by chimeric read.
        if( std::abs((*nodeIter).first - (*breakNode).first) > params->distance){
            continue;
        }
        // Check if there is an edge from current node
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList->find( (*nodeIter).first );
        if( edgeIter==edgeList->end() ){
            continue;
        }
        
        BlockRead totalRead;
        std::vector<ReadVariant> frontRead = this->getBlockRead((*blockVecIter), totalRead, params->crossSNP);
        std::vector<ReadVariant> backRead  = this->getBlockRead((*nextBlockIter), totalRead, params->crossSNP);
        
        std::map<std::string,int> frontReadHP;
        std::map<std::string,int> backReadHP;
        
        // frontReadHP
        // iter all read, determine the haplotype of the read and mark abnormal SNPs
        for(std::vector<ReadVariant>::iterator readIter = frontRead.begin() ; readIter != frontRead.end() ; readIter++ ){
            double hp1Count = 0;
            double hp2Count = 0;
            // loop all variant 
            for( auto variant : (*readIter).variantVec ){
                // get the SNP of the front block
                if( variant.position >= (*breakNode).first ){
                    continue;
                }
                
                PosAllele node = std::make_pair( variant.position , variant.allele+1);
                auto hpIter = subNodeHP->find(node);
                
                if( hpIter != subNodeHP->end() ){
                    if((*subNodeHP)[node]==0){
                        hp1Count++;
                    }
                    else if((*subNodeHP)[node]==1){
                        hp2Count++;
                    }
                } 
                
            }
            // the result of tagging is determined by phased alleles
            if( std::max(hp1Count,hp2Count)/(hp1Count+hp2Count) >= 0.6 && hp1Count+hp2Count >= 1 ){
                int belongHP = ( hp1Count > hp2Count ? 0 : 1 );
                frontReadHP[(*readIter).read_name] = belongHP;
            }
            else{
                frontReadHP[(*readIter).read_name] = -1;
            }
        }    
        // backReadHP
        // iter all read, determine the haplotype of the read and mark abnormal SNPs
        for(std::vector<ReadVariant>::iterator readIter = backRead.begin() ; readIter != backRead.end() ; readIter++ ){
            double hp1Count = 0;
            double hp2Count = 0;
            // loop all variant 
            for( auto variant : (*readIter).variantVec ){
                // get the SNP of the back block
                if( variant.position <= (*nodeIter).first ){
                    continue;
                }
                
                PosAllele node = std::make_pair( variant.position , variant.allele+1);
                std::map<PosAllele,int>::iterator nodePS = bkResult->find(node);
                if( nodePS != bkResult->end() ){
                        if((*subNodeHP)[node]==0)hp1Count++;
                        else hp2Count++;
                }
            }
            
            // the result of tagging is determined by phased alleles
            if( std::max(hp1Count,hp2Count)/(hp1Count+hp2Count) >= 0.6 && hp1Count+hp2Count >= 1 ){
                int belongHP = ( hp1Count > hp2Count ? 0 : 1 );
                backReadHP[(*readIter).read_name] = belongHP;
            }
            else{
                backReadHP[(*readIter).read_name] = -1;
            }
        }   
        
        int rr=0;
        int ra=0;
        int ar=0;
        int aa=0;
        int selectConnect = -1;
        
        // loop all block read, judge each read connect result
        for( std::map<std::string,int>::iterator readIter = totalRead.readVec.begin() ; readIter != totalRead.readVec.end() ; readIter++ ){
            auto frontReadIter = frontReadHP.find((*readIter).first);
            auto backReadIter = backReadHP.find((*readIter).first);
            // this read exists in two blocks
            if( frontReadIter != frontReadHP.end() && backReadIter != backReadHP.end() ){
                if( frontReadHP[(*readIter).first] == 0 && backReadHP[(*readIter).first] == 0 ){
                    rr++;
                }
                else if( frontReadHP[(*readIter).first] == 0 && backReadHP[(*readIter).first] == 1 ){
                    ra++;
                }
                else if( frontReadHP[(*readIter).first] == 1 && backReadHP[(*readIter).first] == 0 ){
                    ar++;
                } 
                else if( frontReadHP[(*readIter).first] == 1 && backReadHP[(*readIter).first] == 1 ){
                    aa++;
                }
                
            }
        }
        
        if( rr + aa > ra + ar ){
            selectConnect = 1;
        }
        else if(rr + aa < ra + ar ){
            selectConnect = 2;
        }
        if ( selectConnect == 1 ){
            PosAllele refStart = std::make_pair((*breakNode).first, 1);
            PosAllele altStart = std::make_pair((*breakNode).first, 2);  
            PosAllele refEnd   = std::make_pair((*nodeIter).first,  1);
            PosAllele altEnd   = std::make_pair((*nodeIter).first,  2);  
            (*bestEdgeConnect)[refEnd] = refStart;
            (*bestEdgeConnect)[altEnd] = altStart;
            connect = true;
        }
        else if ( selectConnect == 2 ){
            PosAllele refStart = std::make_pair((*breakNode).first, 2);
            PosAllele altStart = std::make_pair((*breakNode).first, 1);  
            PosAllele refEnd   = std::make_pair((*nodeIter).first,  1);
            PosAllele altEnd   = std::make_pair((*nodeIter).first,  2);
            (*bestEdgeConnect)[refEnd] = refStart;
            (*bestEdgeConnect)[altEnd] = altStart;
            connect = true;
        }
    }
    return connect;
}
*/
void VairiantGraph::modBridging(){
    
    for(std::map<int,std::vector<int> >::iterator blockVecIter = blockVec->begin() ; blockVecIter != blockVec->end() ; blockVecIter++ ){
        // Number of SNPs checked
        //int sampleNum = 15;
        // get next block
        std::map<int,std::vector<int> >::iterator nextBlockIter = std::next(blockVecIter, 1);
        if( nextBlockIter == blockVec->end() )
            break;
        if( nextBlockIter == blockVecIter )
            break;
        // sort previous block
        std::sort((*blockVecIter).second.rbegin(), (*blockVecIter).second.rend());
        // sort following block
        std::sort((*nextBlockIter).second.begin(), (*nextBlockIter).second.end());
        // The last SNP of the previous block
        std::map<int,ReadBaseMap*>::iterator breakNode = nodeInfo->find((*blockVecIter).second[0]);
        // The last SNP of the following block
        std::map<int,ReadBaseMap*>::iterator NextNodeIter = nodeInfo->find((*nextBlockIter).second[0]);
        // The distance between the two blocks is abnormal, which may be affected by chimeric read.
        if( std::abs((*NextNodeIter).first - breakNode->first) > params->distance){
            continue;
        }
        
        std::map<int,ReadBaseMap*>::iterator fModIter = forwardModNode->begin();
        std::map<int,ReadBaseMap*>::iterator rModIter = reverseModNode->begin();
        
        while( fModIter != forwardModNode->end() && (*fModIter).first < breakNode->first )
            fModIter++;
        
        while( rModIter != reverseModNode->end() && (*rModIter).first < breakNode->first )
            rModIter++;

        int fSelect = subBridge(breakNode->first, NextNodeIter->first, (*forwardModNode));
        int rSelect = subBridge(breakNode->first, NextNodeIter->first, (*reverseModNode));

        if( (fSelect == 1 && rSelect == 2) || (fSelect == 2 && rSelect == 1) ){
            fSelect = -1;
            rSelect = -1;
        }
        
        if( fSelect != -1 && fSelect != -1 ){
            std::cout<< "connect\t" << breakNode->first << "\t~\t" << NextNodeIter->first << "\tdistance\t" << NextNodeIter->first - breakNode->first << "\n\n";
        }

        if ( fSelect == 1 || rSelect == 1 ){
            PosAllele refStart = std::make_pair(breakNode->first, 1);
            PosAllele altStart = std::make_pair(breakNode->first, 2);  
            PosAllele refEnd   = std::make_pair(NextNodeIter->first,  1);
            PosAllele altEnd   = std::make_pair(NextNodeIter->first,  2);  
            (*bestEdgeConnect)[refEnd] = refStart;
            (*bestEdgeConnect)[altEnd] = altStart;
        }
        if ( fSelect == 2 || rSelect == 2 ){
            PosAllele refStart = std::make_pair(breakNode->first, 2);
            PosAllele altStart = std::make_pair(breakNode->first, 1);  
            PosAllele refEnd   = std::make_pair(NextNodeIter->first,  1);
            PosAllele altEnd   = std::make_pair(NextNodeIter->first,  2);
            (*bestEdgeConnect)[refEnd] = refStart;
            (*bestEdgeConnect)[altEnd] = altStart;
        }
        
    }
    
}

int VairiantGraph::subBridge(int breakStart, int breakEnd, std::map<int,ReadBaseMap*> &modNode){

    // store the sub nodes that will be used
    std::map<int,ReadBaseMap*> *subNodeInfo = new std::map<int,ReadBaseMap*>;
    // the list of edges for this subgraph
    std::map<int,VariantEdge*> *subEdgeList = new std::map<int,VariantEdge*>;
    // current snp, haplotype (1 or 2), support snp
    std::map<int, std::map<int,std::vector<int> > > *hpCountMap = new std::map<int, std::map<int,std::vector<int> > >;
    // current snp, result haplotype (1 or 2)
    std::map<int,int> *hpResult = new std::map<int,int>;
    // < block start, <snp in this block> >
    std::map<int,int > *variantBlock = new std::map<int,int>;

    int blockStart = -1;
    int prevSNP = -1;
    int currPos = -1;
    int nextPos = -1;
    int lastConnectPos = -1;

    std::map<int,ReadBaseMap*>::iterator breakStartIter = nodeInfo->find(breakStart);
    std::map<int,ReadBaseMap*>::iterator breakEndIter   = nodeInfo->find(breakEnd);
    
    int prevRange = 4;
    while( breakStartIter != nodeInfo->begin() && prevRange > 0 ){
        breakStartIter = std::prev(breakStartIter, 1);
        prevRange--;
    }
    
    int nextRange = 4;
    while( breakEndIter != nodeInfo->end() && nextRange > 0 ){
        breakEndIter = std::next(breakEndIter, 1);
        nextRange--;
    }


    // build edge process
    // iter all read
    for(std::vector<ReadVariant>::iterator readIter = readVariant->begin() ; readIter != readVariant->end() ; readIter++ ){
        // check start position over than nextPos
        auto lastVariant = (*readIter).variantVec.back();
        if( lastVariant.position < breakStartIter->first ){
            continue;
        }
        // check start position less than breakEnd
        if( (*readIter).reference_start > breakEnd ){
            break;
        }
        
        // --- construct subNodeInfo map ---
        // use for build edge
        ReadVariant tmpRead;
        // iter all variant
        for( auto variant : (*readIter).variantVec ){
            // breakStart and breakEnd will also push in tmpRead
            // check variant position over than breakEnd
            if( variant.position < breakStart ){
                continue;
            }
            // check variant position less than breakEnd
            if( variant.position > breakEnd ){
                break;
            }
            // check the sub nodes that will be used
            // pass break SNP, next SNP and modification position
            auto modIter = modNode.find(variant.position);
            if( variant.position == breakStart || variant.position == breakEnd || modIter != modNode.end() ){
                tmpRead.variantVec.push_back(variant);
            }
            // this position first appear
            auto nodeIter = subNodeInfo->find(variant.position);
            if( nodeIter == subNodeInfo->end() ){
                (*subNodeInfo)[variant.position] = new ReadBaseMap();
            }
            // store the sub nodes
            (*(*subNodeInfo)[variant.position])[(*readIter).read_name] = variant.quality;
        }

        // --- construct subEdgeList map ---
        // iter all pair of snp and construct initial graph
        std::vector<Variant>::iterator variant1Iter = tmpRead.variantVec.begin();
        std::vector<Variant>::iterator variant2Iter = std::next(variant1Iter,1);
        while(variant1Iter != tmpRead.variantVec.end() && variant2Iter != tmpRead.variantVec.end() ){
            // create new edge if not exist
            std::map<int,VariantEdge*>::iterator posIter = subEdgeList->find((*variant1Iter).position);
            if( posIter == subEdgeList->end() )
                (*subEdgeList)[(*variant1Iter).position] = new VariantEdge((*variant1Iter).position);
            // add edge process
            for(int nextNode = 0 ; nextNode < params->connectAdjacent ; nextNode++){
                // prevent SNP connect SNP
                if( (*variant1Iter).position == breakStart && (*variant2Iter).position == breakEnd )
                    continue;
                // this position is mod
                if( (*variant1Iter).allele == 0 )
                    (*subEdgeList)[(*variant1Iter).position]->ref->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).read_name);
                // this position is nonMod
                if( (*variant1Iter).allele == 1 )
                    (*subEdgeList)[(*variant1Iter).position]->alt->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).read_name);
                
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
    
    //---------------------------------------------------------
    
    for(auto posIter = subNodeInfo->begin() ; posIter != subNodeInfo->end() ; posIter++ ){
        std::cout<< (*posIter).first << " ";
    }
    
    std::cout<<"\n";  
    
    for(std::map<int,VariantEdge*>::iterator edgeIter = subEdgeList->begin(); edgeIter != subEdgeList->end() ; edgeIter++){
        std::vector<std::string> re1 = (*edgeIter).second->ref->showEdge("ref_"+std::to_string((*edgeIter).first));
        std::vector<std::string> re2 = (*edgeIter).second->alt->showEdge("alt_"+std::to_string((*edgeIter).first));

        for( auto edge : re1){
            std::cout<< edge << "\n";
        }
        for( auto edge : re2){
            std::cout<< edge << "\n";
        }
    }
    std::cout<<"\n";
    //---------------------------------------------------------
    
    
    // try to filter noise position
    
    bool skipEnd = false;
    // pos, < allele, read count >
    std::map<int,std::map<int,int>> *innerEdgeCount = new std::map<int,std::map<int,int>>;
    
    for(std::map<int,VariantEdge*>::iterator edgeIter = subEdgeList->begin(); edgeIter != subEdgeList->end() ; edgeIter++){
        // < position, allele >
        // get all positions connected from the current position
        std::vector<std::pair<int,int>> refConnetPos = (*edgeIter).second->ref->getConnectPos();
        std::vector<std::pair<int,int>> altConnetPos = (*edgeIter).second->alt->getConnectPos();
        // count each position's inner edge
        // the position connect from ref allele
        for(auto posIter : refConnetPos )
            (*innerEdgeCount)[posIter.first][posIter.second]++;
        // the position connect from alt allele
        for(auto posIter : altConnetPos )
            (*innerEdgeCount)[posIter.first][posIter.second]++;
    }
    
    for(auto innerEdge = innerEdgeCount->begin() ; innerEdge != innerEdgeCount->end() ; innerEdge++ ){
        for( auto posIter = innerEdge->second.begin() ; posIter != innerEdge->second.end() ; posIter++ ){
            std::cout <<"count\t"<< innerEdge->first << "\t" << posIter->first << "\t" << posIter->second << "\n";
            // because default connectAdjacent is 6
            // if this position has more than 10 inner edges,
            // it means both reference and alternative connect to this position 
            if( posIter->second >= 10 ){
                std::cout<<"erase\t" << innerEdge->first <<"\n";
                
                subNodeInfo->erase( innerEdge->first );

                if( innerEdge->first == breakEnd ){
                    skipEnd = true;
                }
            }
        }
    }
    
    delete innerEdgeCount;
    
    for(std::map<int,VariantEdge*>::iterator edgeIter = subEdgeList->begin(); edgeIter != subEdgeList->end() ; edgeIter++){
		(*edgeIter).second->ref->deletEdge((*edgeIter).first);
		(*edgeIter).second->alt->deletEdge((*edgeIter).first);
    }
    

    // run graph algorithm 
    if(!skipEnd)
    {
        // Visit all position and assign SNPs to haplotype
        // Only one of the two alleles needs to be used for each SNP
        for(std::map<int,ReadBaseMap*>::iterator nodeIter = subNodeInfo->begin() ; nodeIter != subNodeInfo->end() ; nodeIter++ ){
            
            // check next position
            std::map<int,ReadBaseMap*>::iterator nextNodeIter = std::next(nodeIter, 1);
            if( nextNodeIter == subNodeInfo->end() ){
                int h1 = (*hpCountMap)[nodeIter->first][1].size();
                int h2 = (*hpCountMap)[nodeIter->first][2].size();
                
                if( h1 == 0 && h2 == 0 ){
                    (*variantBlock)[nodeIter->first] = nodeIter->first;
                    (*hpResult)[nodeIter->first] = 1;
                }
                else if( h1 > h2 || h1 < h2 ){
                    int currHP = ( h1 > h2 ? 1 : 2 );
                    (*hpResult)[nodeIter->first] = currHP;
                    (*variantBlock)[nodeIter->first] = blockStart;
                }
                break;
            }
            
            prevSNP = currPos;
            currPos = nodeIter->first;
            nextPos = nextNodeIter->first;
            
            // check distance between two position 
            if(std::abs(nextPos-currPos) > params->distance ){
                continue;
            }
            
            // get the number of HP1 and HP2 supported reference allele
            int h1 = (*hpCountMap)[currPos][1].size();
            int h2 = (*hpCountMap)[currPos][2].size();
            std::cout<< "h1:h2\t" << currPos << "\t" << h1 << "\t" << h2 << "\n";
            double phasedConfidence   = (double)std::max(h1,h2)/(double)(h1+h2);    
            
            // new block, set this position as block start 
            if( h1 == 0 && h2 == 0 ){
                
                // No new blocks should be created if the next SNP has already been picked up
                if( currPos < lastConnectPos ){
                    continue;
                }
                
                blockStart=currPos;
                (*variantBlock)[currPos]=blockStart;
                (*hpResult)[currPos] = 1;
            }
            // ignore poorly classified SNPs
            else if( 1 - params->confidentHaplotype <= phasedConfidence  && phasedConfidence <= params->confidentHaplotype ){
                // check previos SNP in edge list
                std::map<int,VariantEdge*>::iterator edgeIter = subEdgeList->find( prevSNP );
                if( edgeIter != subEdgeList->end() ){
                    // record none haplotype
                    std::pair<PosAllele,PosAllele> preResult = (*edgeIter).second->findBestEdgePair(currPos, params->isONT, params->readsThreshold, true, false);  
                    if( preResult.first.second == 1 || preResult.first.second == 2 ){
                        (*hpResult)[currPos] = 0;
                        (*variantBlock)[currPos] = blockStart;
                    }
                }
            }
            else{
                if( h1 > h2 || h1 < h2 ){
                    int currHP = ( h1 > h2 ? 1 : 2 );
                    (*hpResult)[currPos] = currHP;
                    (*variantBlock)[currPos] = blockStart;
                }
            }
            
            // Check if there is no edge from current node
            std::map<int,VariantEdge*>::iterator edgeIter = subEdgeList->find( currPos );
            if( edgeIter==subEdgeList->end() )
                continue;
            // check connect between surrent SNP and next n SNPs
            for(int i = 0 ; i < params->connectAdjacent ; i++ ){
                // consider reads from the currnt SNP and the next (i+1)'s SNP
                std::pair<PosAllele,PosAllele> tmp = edgeIter->second->findBestEdgePair(nextNodeIter->first, params->isONT, 1, false, false);
                
                //
                // -1 : no connect  
                //  1 : the haplotype of next (i+1)'s SNP are same as previous
                //  2 : the haplotype of next (i+1)'s SNP are different as previous
                if( tmp.first.second != -1 ){
                    // record the haplotype resut of next (i+1)'s SNP
                    if( (*hpResult)[currPos] == 1 ){
                        if( tmp.first.second == 1 ){
                            (*hpCountMap)[nextNodeIter->first][1].push_back(currPos);
                            //std::cout<<"connect 11\t" << currPos << " -> " <<  nextNodeIter->first << "\n";
                        }
                        if( tmp.first.second == 2 ){
                            (*hpCountMap)[nextNodeIter->first][2].push_back(currPos);
                            //std::cout<<"connect 12\t" << currPos << " -> " <<  nextNodeIter->first << "\n";
                        }
                    }
                    if( (*hpResult)[currPos]==2 ){
                        if( tmp.first.second == 1 ){
                            (*hpCountMap)[nextNodeIter->first][2].push_back(currPos);
                            //std::cout<<"connect 21\t" << currPos << " -> " <<  nextNodeIter->first << "\n";
                        }
                        if( tmp.first.second == 2 ){
                            (*hpCountMap)[nextNodeIter->first][1].push_back(currPos);
                            //std::cout<<"connect 22\t" << currPos << " -> " <<  nextNodeIter->first << "\n";
                        }
                    }
                    
                    
                    lastConnectPos = nextNodeIter->first;
                }
                
                std::cout<<"status\t" << currPos << " -> " <<  nextNodeIter->first 
                             << "\t" << tmp.first.second << "\t" << (*hpResult)[currPos] << "\t|\t"<< (*hpCountMap)[nextNodeIter->first][1].size() << "\t" <<(*hpCountMap)[nextNodeIter->first][2].size() << "\n";
                    
                nextNodeIter++;
                if( nextNodeIter == subNodeInfo->end() )
                    break;
            }
        }
    }
/*         
    
    for(auto posIter = subNodeInfo->begin() ; posIter != subNodeInfo->end() ; posIter++ ){
        int h1 = (*hpCountMap)[(*posIter).first ][1].size();
        int h2 = (*hpCountMap)[(*posIter).first][2].size();
        std::cout<<"HP\t"<< (*posIter).first << "\tHP:" << h1 << "\t" << h2 << "\tblock: " << (*variantBlock)[(*posIter).first] << "\n";
        
    }
    
*/
    
    int startPS = (*variantBlock)[breakStart];
    int endPS   = (*variantBlock)[breakEnd];
    int startHP = (*hpResult)[breakStart];
    int endHP   = (*hpResult)[breakEnd];
    int select = -1;

    if( startPS == endPS && !skipEnd ){
        if( startHP == endHP){
            select = 1;
        }
        else{
            select = 2;
        }
    }

    std::cout<< "containMethyl\t" << "\t" << select << "\t" << subNodeInfo->size() - 2 << "\n";

    
    //std::cout<< "result:\t" << select << "\t"<< breakStart << "," << startPS << ","<< startHP << "\t" 
    //                        << breakEnd << "," << endPS << ","<< endHP << "\n";
    
    ///---------------------------------------------------------
    for( auto edgeIter = subEdgeList->begin() ; edgeIter != subEdgeList->end() ; edgeIter++ ){
        edgeIter->second->ref->destroy();
        edgeIter->second->alt->destroy();
        delete edgeIter->second->ref;
        delete edgeIter->second->alt;
    }
    
    for( auto nodeIter = subNodeInfo->begin() ; nodeIter != subNodeInfo->end() ; nodeIter++ ){
        delete nodeIter->second;
    }
    
    delete subEdgeList;
    delete subNodeInfo;
    delete hpCountMap;
    delete hpResult;
    delete variantBlock;

    return select;
}

void VairiantGraph::storeResultPath(){
    // Initialize the parameters of the stored results
    posAppear->clear();
    blockStart->clear();
    bkResult->clear();
    subNodeHP->clear();
    blockVec->clear();

    // loop all best edge, try to restore the full connect path
    for(std::map<PosAllele,PosAllele>::iterator edgeIter = bestEdgeConnect->begin() ; edgeIter != bestEdgeConnect->end() ; edgeIter++ ){
        PosAllele A = (*edgeIter).first;
        PosAllele B = (*edgeIter).second;
        
        PosAllele aOtherSide = std::make_pair(A.first, (A.second != 1 ? 1 : 2));
        PosAllele bOtherSide = std::make_pair(B.first, (B.second != 1 ? 1 : 2));

        // in order to connect block path
        std::map<PosAllele,int>::iterator aBlockIter = bkResult->find(A);
        std::map<PosAllele,int>::iterator bBlockIter = bkResult->find(B);

        std::map<int,int>::iterator aAppearIter = posAppear->find(A.first);
        std::map<int,int>::iterator bAppearIter = posAppear->find(B.first);

        std::map<int,int>::iterator aStartIter = blockStart->find(A.first);
        std::map<int,int>::iterator bStartIter = blockStart->find(B.first);
        
        int minBlockStart = std::min(A.first+1, B.first+1);

        if( aStartIter != blockStart->end() )
            minBlockStart = std::min(minBlockStart, (*aStartIter).second );
        if( bStartIter != blockStart->end() )
            minBlockStart = std::min(minBlockStart, (*bStartIter).second );
        
        // two node doesn't appear, block start
        if( aBlockIter == bkResult->end() && bBlockIter == bkResult->end() ){
            (*bkResult)[A] = minBlockStart;
            (*bkResult)[B] = minBlockStart;
            // add HP tag
            if( aAppearIter == posAppear->end() && bAppearIter == posAppear->end() ){
                (*subNodeHP)[A] = 0;
                (*subNodeHP)[B] = 0;
                (*subNodeHP)[aOtherSide] = 1;
                (*subNodeHP)[bOtherSide] = 1;
            }
        }
        // node A doesn't appear
        else if( aBlockIter == bkResult->end() && bBlockIter != bkResult->end() ){
            (*bkResult)[A] = std::min(minBlockStart, (*bkResult)[B]);
            (*bkResult)[aOtherSide] = std::min(minBlockStart, (*bkResult)[B]);
            
            (*subNodeHP)[A] = (*subNodeHP)[B];
            (*subNodeHP)[aOtherSide] = ( (*subNodeHP)[A] != 0 ? 0 : 1);//subNodeHP[bOtherSide];
        }
        // node B doesn't appear
        else if( aBlockIter != bkResult->end() && bBlockIter == bkResult->end() ){
            (*bkResult)[B] = std::min(minBlockStart, (*bkResult)[A]);
            (*bkResult)[bOtherSide] = std::min(minBlockStart, (*bkResult)[A]);
            
            (*subNodeHP)[B] = (*subNodeHP)[A];
            (*subNodeHP)[bOtherSide] = ( (*subNodeHP)[B] != 0 ? 0 : 1);//subNodeHP[aOtherSide];
        }
        else{
            // For diploid,
            // In most cases, there will be two sub-edges between every two positions.
            // And the phasing of two positions has been completed on the first sub-edge.
        }

        (*posAppear)[A.first] = 1;
        (*posAppear)[B.first] = 1;
        
        (*blockStart)[A.first] = minBlockStart;
        (*blockStart)[B.first] = minBlockStart;
        
        if( aAppearIter == posAppear->end())
            (*blockVec)[minBlockStart].push_back(A.first);
        if( bAppearIter == posAppear->end())
            (*blockVec)[minBlockStart].push_back(B.first);
    }
}

VairiantGraph::VairiantGraph(std::string &in_ref, PhasingParameters &in_params){
    params=&in_params;
    ref=&in_ref;
    
    nodeInfo = new std::map<int,ReadBaseMap*>;
    edgeList = new std::map<int,VariantEdge*>;
    bestEdgeConnect = new std::map<PosAllele,PosAllele>;
    posAppear = new std::map<int,int>;
    blockStart = new std::map<int,int>;
    forwardModNode = new std::map<int,ReadBaseMap*>;
    reverseModNode = new std::map<int,ReadBaseMap*>;
    bkResult = new std::map<PosAllele,int>;
    subNodeHP = new std::map<PosAllele,int>;
    blockVec = new std::map<int,std::vector<int>> ;
    svPosition = new std::map<int,int>;
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
    delete bestEdgeConnect;
    delete posAppear;
    delete blockStart;
    delete forwardModNode;
    delete reverseModNode;
    delete bkResult;
    delete subNodeHP;
    delete blockVec;  
    delete readHpMap;
}

void VairiantGraph::addEdge(std::vector<ReadVariant> &in_readVariant){
    readVariant = &in_readVariant;
    // iter all read
    for(std::vector<ReadVariant>::iterator readIter = in_readVariant.begin() ; readIter != in_readVariant.end() ; readIter++ ){
        ReadVariant tmpRead;
        //std::cout<< "read\t" << (*readIter).source_id << "\t" << (*readIter).read_name << "\t" << (*readIter).reference_start << "\n";
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
            
            // assign SV quality
            if( variant.quality == -1 ){
                (*svPosition)[variant.position] = 1;
                if( variant.allele == 1 ){
                    // SV calling
                    variant.quality = 60; 
                }
                else{
                    // assume ref
                    variant.quality = 30;
                }
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
                    (*edgeList)[(*variant1Iter).position]->ref->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).read_name);
                // this allele support alt
                if( (*variant1Iter).allele == 1 )
                    (*edgeList)[(*variant1Iter).position]->alt->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).read_name);
                
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
            PosAllele node = std::make_pair( variant.position , variant.allele+1);
            std::map<PosAllele,int>::iterator nodePS = bkResult->find(node);

            if( nodePS != bkResult->end() ){
                if( (*bkResult)[node] != 0 ){
                    if((*subNodeHP)[node]==0)refCount++;
                    else altCount++;
                }
            }
        }

        // tag high confident reads
        if( refCount == altCount ){
            (*readHpMap)[(*readIter).read_name] = -1;
        }
        else if( std::max(refCount,altCount)/(refCount+altCount) >= params->readConfidence ){
            // tag read with the corresponding haplotype
            int belongHP = ( refCount > altCount ? 0 : 1 );
            (*readHpMap)[(*readIter).read_name] = belongHP;
        }
        else{
            (*readHpMap)[(*readIter).read_name] = -1;
        }
    }

    int totalBase = 0;
    std::map<int,int> *coverBase = new std::map<int,int>;
    
    // count each haplotype's allele on SNP locus
    for(auto readIter = (*readVariant).begin() ; readIter != (*readVariant).end() ; readIter++ ){
        // Skip reads with insufficient confidence
        auto hpIter = readHpMap->find((*readIter).read_name);
        if( hpIter == readHpMap->end() || (*readHpMap)[(*readIter).read_name] == -1 ){
            continue;
        }
        
        int readHP = (*readHpMap)[(*readIter).read_name];
        for(auto variantIter = (*readIter).variantVec.begin() ; variantIter != (*readIter).variantVec.end() ; variantIter++ ){
            if( ( (*variantIter).allele == 0 || (*variantIter).allele == 1 ) && (*variantIter).quality >= 0 ){
                int allele = (*variantIter).allele;
                (*hpAlleleCountMap)[readHP][(*variantIter).position][allele]++;
                
                totalBase++;
                (*coverBase)[(*variantIter).position] = 1;
            }
        }
    }

    if( totalBase != 0 && coverBase->size() != 0 ){
        int avgDepth = totalBase/coverBase->size();
        double snpConfidenceThreshold = 0.6;
        // when the read coverage is low, the higher threshold can reduce the impact of few incorrect reads.
        if( avgDepth < 20 ){
            snpConfidenceThreshold = 0.8;
        }
        else if( avgDepth < 30 ){
            snpConfidenceThreshold = 0.7;
        }
        

        std::cerr<< "AvgDepth " << avgDepth << " ... ";

        subNodeHP->clear();
        
        std::map<int,std::map<int,int>> hpAllele;
        // reassign allele result
        for(auto nodeIter = nodeInfo->begin() ; nodeIter != nodeInfo->end() ; nodeIter++ ){
            int position = nodeIter->first;
            PosAllele A = std::make_pair(position, 1);
            PosAllele aOtherSide = std::make_pair(position, 2);
            
            double hp1Ref = (*hpAlleleCountMap)[0][position][0];
            double hp1Alt = (*hpAlleleCountMap)[0][position][1];
            double hp2Ref = (*hpAlleleCountMap)[1][position][0];
            double hp2Alt = (*hpAlleleCountMap)[1][position][1];
            double result1reads = hp1Ref + hp2Alt;
            double result2reads = hp2Ref + hp1Alt;
            double resultConfidence = std::max(result1reads, result2reads) / (result1reads + result2reads);
            
            int hp1Result = -1;
            int hp2Result = -1;
            
            std::cout<<position<<"\t"<<resultConfidence<<"\t"<<snpConfidenceThreshold<<"\t"<<hp1Ref<<"\t"<<hp1Alt<<"\t"<<hp2Ref<<"\t"<<hp2Alt<<"\t"<<result1reads<<"\t"<<result2reads<<"\n";
            
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
                (*subNodeHP)[A] = hp1Result;
                (*subNodeHP)[aOtherSide] = hp2Result;
            }
            else{
                bkResult->erase(A);
                bkResult->erase(aOtherSide);
            }
        }
    }

    delete coverBase;
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
    for( std::map<int,ReadBaseMap*>::iterator nodeIter = nodeInfo->begin() ; nodeIter != nodeInfo->end() ; nodeIter++ ){
        
        PhasingElement tmp;
        
        PosAllele ref = std::make_pair( nodeIter->first , 1);
        PosAllele alt = std::make_pair( nodeIter->first , 2);
        
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
            std::string key = chrName + "_" + std::to_string( nodeIter->first );
            result[key] = tmp;
        }
    }
}

std::map<std::string,int>* VairiantGraph::getReadHP(){
    return readHpMap;
}

int VairiantGraph::totalNode(){
    return nodeInfo->size();
}

void VairiantGraph::connectTest(std::string chrName, std::map<int,int> &passPosition){


    // position, number of clear edge, number of dirty (unclear) edge
    std::map<int, std::pair<int, int> > *innerEdge = new std::map<int, std::pair<int, int> >;
    std::map<int, std::pair<int, int> > *outerEdge = new std::map<int, std::pair<int, int> >;

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

            double majorRatio = (double)std::max(tmp.first,tmp.second)/(double)(tmp.first+tmp.second);

            //if( majorRatio >= 0.8 && tmp.first + tmp.second >= 4){
                        /*if( majorRatio >= 0.8 && tmp.first + tmp.second >= (((*(*nodeIter).second).size() + (*(*nextNodeIter).second).size())/4)){
                (*outerEdge)[currPos].first++;
                (*innerEdge)[nextPos].first++;
            }
            else{
                (*outerEdge)[currPos].second++;
                (*innerEdge)[nextPos].second++;
            }*/

            std::cout<< chrName << "\t" <<"RRxRA\t" << currPos << "\t->\t" << nextPos << "\t | RR\t"<< tmp.first << "\tRA\t"<< tmp.second << "\t|\t" << majorRatio << "\t" << (((*(*nodeIter).second).size() + (*(*nextNodeIter).second).size())/3) << "\n";


            //if( majorRatio >= 0.95 && tmp.first + tmp.second > (((*(*nodeIter).second).size() + (*(*nextNodeIter).second).size())/4)){
            if( majorRatio >= 0.9 &&  tmp.first + tmp.second > (((*(*nodeIter).second).size() + (*(*nextNodeIter).second).size())/4) && tmp.first + tmp.second > 6 ){
                    passPosition[currPos] = 1;
                    passPosition[nextPos] = 1;
            }
            nextNodeIter++;
            if( nextNodeIter == nodeInfo->end() )
                break;
        }
    }

    // check all variant inner and outer edge

    /*for(std::map<int,ReadBaseMap*>::iterator nodeIter = nodeInfo->begin() ; nodeIter != nodeInfo->end() ; nodeIter++ ){
    int currPos = nodeIter->first;

    double clearEdge = (*innerEdge)[currPos].first  + (*outerEdge)[currPos].first;
    double dirtyEdge = (*innerEdge)[currPos].second + (*outerEdge)[currPos].second;
    double clearRatio = clearEdge/(double)(clearEdge+dirtyEdge);

    std::cout<< "edge\t" << currPos << "\t| innerEdge\t" << (*innerEdge)[currPos].first << "\t" << (*innerEdge)[currPos].second
                        << "\t| outerEdge\t" << (*outerEdge)[currPos].first << "\t" << (*outerEdge)[currPos].second
                        << "\t" << clearRatio << "\n";

            if( clearRatio > 0.25){
                    passPosition[currPos] = 1;
            }
    }*/

    delete innerEdge;
    delete outerEdge;

}

void VairiantGraph::phasingProcess(){

    this->edgeConnectResult();
    this->storeResultPath();
    
    this->modBridging();
    this->storeResultPath();
    
    /*
    while(true){
        bool connect = this->blockPhaseing();
        this->storeResultPath();
        if(!connect){
            break;
        }
    }
    */
    this->readCorrection();  
}