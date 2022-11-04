#include "PhasingGraph.h"

//SubEdge

SubEdge::SubEdge():readCount(0){ 
}

void SubEdge::addSubEdge(int currentQuality, Variant connectNode, std::string readName){
    if(connectNode.allele == 0 ){
        refRead[connectNode.position].push_back(readName);
        // quality sum
        std::map<int, int>::iterator rqIter = refQuality.find(connectNode.position);
        if( rqIter == refQuality.end() ){
            refQuality[connectNode.position] = currentQuality + connectNode.quality;
            refReadCount[connectNode.position] = 1;
        }
        else{
            refQuality[connectNode.position] += currentQuality + connectNode.quality;
            refReadCount[connectNode.position]++;
        }
    }
    if(connectNode.allele == 1 ){
        altRead[connectNode.position].push_back(readName);
        // quality sum
        std::map<int, int>::iterator aqIter = altQuality.find(connectNode.position);
        if( aqIter == altQuality.end() ){
            altQuality[connectNode.position] = currentQuality + connectNode.quality;
            altReadCount[connectNode.position] = 1;
        }
        else{
            altQuality[connectNode.position] += currentQuality + connectNode.quality;
            altReadCount[connectNode.position]++;
        }
    }
    readCount++;
}

std::vector<std::string> SubEdge::showEdge(std::string message){
    std::vector<std::string> result;
    for(std::map<int, std::vector<std::string> >::iterator edgeIter = refRead.begin() ; edgeIter != refRead.end() ; edgeIter++ ){
        result.push_back(message +" -> ref_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second.size()) + "];");
    }
    for(std::map<int, std::vector<std::string> >::iterator edgeIter = altRead.begin() ; edgeIter != altRead.end() ; edgeIter++ ){
        result.push_back(message +" -> alt_" + std::to_string((*edgeIter).first) + "[label=" + std::to_string((*edgeIter).second.size()) + "];");
    }
    return result;
}

std::pair<int,int> SubEdge::BestPair(int targetPos){
    // default, no read support target position
    int refResult = 0;
    int altResult = 0;
    // find the edge between local and target position 
    std::map<int, std::vector<std::string> >::iterator refEdgeIter = refRead.find(targetPos);
    std::map<int, std::vector<std::string> >::iterator altEdgeIter = altRead.find(targetPos);
    // local connect target REF allele
    if( refEdgeIter != refRead.end() ){
        refResult = (*refEdgeIter).second.size();
    }
    // local connect target ALT allele
    if( altEdgeIter != altRead.end() ){
        altResult = (*altEdgeIter).second.size();   
    }            
    return std::make_pair( refResult, altResult );
}

void SubEdge::getReadVariant(std::map< std::string,std::map<int,int> > &readVariantMap){

    for(std::map<int, std::vector<std::string> >::iterator edgeIter = refRead.begin() ; edgeIter != refRead.end() ; edgeIter++ ){
        for(std::vector<std::string>::iterator readIter = (*edgeIter).second.begin() ; readIter != (*edgeIter).second.end() ; readIter++ ){
            readVariantMap[(*readIter)][(*edgeIter).first]=0;
        }
    }
    
    for(std::map<int, std::vector<std::string> >::iterator edgeIter = altRead.begin() ; edgeIter != altRead.end() ; edgeIter++ ){
        for(std::vector<std::string>::iterator readIter = (*edgeIter).second.begin() ; readIter != (*edgeIter).second.end() ; readIter++ ){
            readVariantMap[(*readIter)][(*edgeIter).first]=1;
        }
    }
}

int SubEdge::getQuality(PosAllele targetPos){
    // target is Ref allele
    if( targetPos.second == 1 ){
        std::map<int, int>::iterator qIter = refQuality.find(targetPos.first);
        if( qIter == refQuality.end() )
            return 0;
        else
            return refQuality[targetPos.first];
    }
    // target is Alt allele
    if( targetPos.second == 2 ){
        std::map<int, int>::iterator qIter = altQuality.find(targetPos.first);
        if( qIter == altQuality.end() )
            return 0;
        else
            return altQuality[targetPos.first];
    }
    return 0;
}

int SubEdge::getAvgQuality(PosAllele targetPos){
    // target is Ref allele
    if( targetPos.second == 1 ){
        std::map<int, int>::iterator qIter = refQuality.find(targetPos.first);
        if( qIter == refQuality.end() )
            return 0;
        else
            return refQuality[targetPos.first]/refReadCount[targetPos.first];
    }
    // target is Alt allele
    if( targetPos.second == 2 ){
        std::map<int, int>::iterator qIter = altQuality.find(targetPos.first);
        if( qIter == altQuality.end() )
            return 0;
        else
            return altQuality[targetPos.first]/altReadCount[targetPos.first];
    }
    return 0;
}

int SubEdge::getRefReadCount(int targetPos){
    std::map<int, int>::iterator posIter = refReadCount.find(targetPos);
    if( posIter != refReadCount.end() ){
        return refReadCount[targetPos];
    }
    return 0;
}

int SubEdge::getAltReadCount(int targetPos){
    std::map<int, int>::iterator posIter = altReadCount.find(targetPos);
    if( posIter != altReadCount.end() ){
        return altReadCount[targetPos];
    }
    return 0;
}

VariantEdge::VariantEdge(int inCurrPos){
    currPos = inCurrPos;
    alt = new SubEdge();
    ref = new SubEdge();
}

//VariantEdge
std::pair<PosAllele,PosAllele> VariantEdge::findBestEdgePair(int targetPos, bool isONT, double diffRatioThreshold, double svDiffRatioThreshold, bool containSV, bool repair){
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
    std::map<int, std::map<int,std::vector<int> > > hpCountMap;
    // current snp, result haplotype (1 or 2)
    std::map<int,int> hpResult;
    // current snp, result haplotype (1 or 2)
    std::map<int,int> inconsistentSnpMap;
    // < block start, <snp in this block> >
    std::map<int,std::vector<int> > phasedBlocks;
    int blockStart = -1;
    
    int currPos = -1;
    int nextPos = -1;
    
    int lastConnectPos = -1;
    
    // Visit all position and assign SNPs to haplotype
    // Only one of the two alleles needs to be used for each SNP
    for(std::map<int,ReadBase>::iterator nodeIter = nodeInfo.begin() ; nodeIter != nodeInfo.end() ; nodeIter++ ){
        // check next position
        std::map<int,ReadBase>::iterator nextNodeIter = std::next(nodeIter, 1);
        if( nextNodeIter == nodeInfo.end() ){
             break;
        }
        
        currPos = (*nodeIter).first;
        nextPos = (*nextNodeIter).first;
        
        // check distance between two position 
        if(std::abs(nextPos-currPos) > params->distance ){
            continue;
        }
        
        // get the number of HP1 and HP2 supported reference allele
        int h1 = hpCountMap[currPos][1].size();
        int h2 = hpCountMap[currPos][2].size();
        
        double haplotypeConsistentRatio   = (double)std::max(h1,h2)/(double)(h1+h2);
        double haplotypeInconsistentRatio = (double)std::min(h1,h2)/(double)(h1+h2);
        
        // check and record inconsistent snp
        // using haplotype ratio to determine whether inconsistency occurs
        if( h1!=0 && h2!=0 && haplotypeInconsistentRatio <= params->judgeInconsistent ){
            // calculate h1 ratio
            double h1ratio = (double)h1/(double)(h1+h2);

            // support h1 reads are inconsistent
            if( h1ratio == haplotypeInconsistentRatio ){
                for( unsigned int i = 0 ; i < hpCountMap[currPos][1].size() ; i++ ){
                    inconsistentSnpMap[hpCountMap[currPos][1][i]]++;
                }
            }
            // support h2 reads are inconsistent
            else{
                for( unsigned int i = 0 ; i < hpCountMap[currPos][2].size() ; i++ ){
                    inconsistentSnpMap[hpCountMap[currPos][2][i]]++;
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
            phasedBlocks[blockStart].push_back(currPos);
            hpResult[currPos] = 1;
        }
        
        // ignore poorly classified SNPs
        // this SNP can not judge by forword n SNPs
        // using the reads between previous SNP to judge current SNP haplotype
        else if( 1 - params->confidentHaplotype <= haplotypeConsistentRatio  && haplotypeConsistentRatio <= params->confidentHaplotype ){
            hpResult[currPos] = 0;
        }
        else{
            if( h1 > h2 || h1 < h2 ){
                int currHP = ( h1 > h2 ? 1 : 2 );
                hpResult[currPos] = currHP;
                phasedBlocks[blockStart].push_back(currPos);
            }
        }
        // Check if there is no edge from current node
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList.find( currPos );
        if( edgeIter==edgeList.end() )
            continue;
        // check connect between surrent SNP and next n SNPs
        for(int i = 0 ; i < params->connectAdjacent ; i++ ){
            // consider reads from the currnt SNP and the next (i+1)'s SNP
            std::pair<PosAllele,PosAllele> tmp = (*edgeIter).second->findBestEdgePair((*nextNodeIter).first, params->isONT, 1, 1, false, false);
            // -1 : no connect  
            //  1 : the haplotype of next (i+1)'s SNP are same as previous
            //  2 : the haplotype of next (i+1)'s SNP are different as previous
            if( tmp.first.second != -1 ){
                // record the haplotype resut of next (i+1)'s SNP
                if( hpResult[currPos] == 1 ){
                    if( tmp.first.second == 1 ){
                        hpCountMap[(*nextNodeIter).first][1].push_back(currPos);
                    }
                    if( tmp.first.second == 2 ){
                        hpCountMap[(*nextNodeIter).first][2].push_back(currPos);
                    }
                }
                if( hpResult[currPos]==2 ){
                    if( tmp.first.second == 1 ){
                        hpCountMap[(*nextNodeIter).first][2].push_back(currPos);
                    }
                    if( tmp.first.second == 2 ){
                        hpCountMap[(*nextNodeIter).first][1].push_back(currPos);
                    }
                }
                if( params->generateDot ){
                    std::string e1 = std::to_string(currPos+1) + ".1\t->\t" + std::to_string(tmp.first.first+1) + "." + std::to_string(tmp.first.second);
                    std::string e2 = std::to_string(currPos+1) + ".2\t->\t" + std::to_string(tmp.second.first+1) + "." + std::to_string(tmp.second.second);

                    dotResult.push_back(e1);
                    dotResult.push_back(e2);
                }
                
                lastConnectPos = (*nextNodeIter).first;
            }
            nextNodeIter++;
            if( nextNodeIter == nodeInfo.end() )
                break;
        }
    }

    // Correction inconsistent SNPs
    std::vector<std::string>::iterator edgeIter = dotResult.begin();
    for(std::map<int,int>::iterator snpIter = inconsistentSnpMap.begin();snpIter != inconsistentSnpMap.end();snpIter++){
        if( (*snpIter).second >= params->inconsistentThreshold ){
            // change haplotype result
            hpResult[(*snpIter).first] = (hpResult[(*snpIter).first] == 1 ? 2 : 1);
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
    for(auto blockIter = phasedBlocks.begin() ; blockIter != phasedBlocks.end() ; blockIter++ ){
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
            
            if( hpResult[(*currIter)] == 0 || hpResult[(*nextIter)] == 0 ){
                
            }
            else if( hpResult[(*currIter)] == hpResult[(*nextIter)] ){
                refEnd = std::make_pair((*nextIter), 1);
                altEnd = std::make_pair((*nextIter), 2);
            }
            else{
                refEnd = std::make_pair((*nextIter), 2);
                altEnd = std::make_pair((*nextIter), 1);
            }
            bestEdgeConnect[refEnd] = refStart;
            bestEdgeConnect[altEnd] = altStart;
        }
    }
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
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList.find( (*posIter) );
        if( edgeIter==edgeList.end() )
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

bool VairiantGraph::blockPhaseing(){
    
    bool connect = false;
    
    // loop all block, prepare to create new edges between two block
    // PS number, < position vector >
    for(std::map<int,std::vector<int> >::iterator blockVecIter = blockVec.begin() ; blockVecIter != blockVec.end() ; blockVecIter++ ){
        // Number of SNPs checked
        //int sampleNum = 15;
        // get next block
        std::map<int,std::vector<int> >::iterator nextBlockIter = std::next(blockVecIter, 1);
        if( nextBlockIter == blockVec.end() )
            break;
        if( nextBlockIter == blockVecIter )
            break;
        // sort previous block
        std::sort((*blockVecIter).second.rbegin(), (*blockVecIter).second.rend());
        // sort following block
        std::sort((*nextBlockIter).second.begin(), (*nextBlockIter).second.end());
        // The last SNP of the previous block
        std::map<int,ReadBase>::iterator breakNode = nodeInfo.find((*blockVecIter).second[0]);
        // The last SNP of the following block
        std::map<int,ReadBase>::iterator nodeIter = nodeInfo.find((*nextBlockIter).second[0]);
        // The distance between the two blocks is abnormal, which may be affected by chimeric read.
        if( std::abs((*nodeIter).first - (*breakNode).first) > params->distance){
            continue;
        }
        // Check if there is an edge from current node
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList.find( (*nodeIter).first );
        if( edgeIter==edgeList.end() ){
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
                
                auto hpIter = subNodeHP.find(node);
                
                if( hpIter != subNodeHP.end() ){
                    if(subNodeHP[node]==0){
                        hp1Count++;
                    }
                    else if(subNodeHP[node]==1){
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
                std::map<PosAllele,int>::iterator nodePS = bkResult.find(node);
                if( nodePS != bkResult.end() ){
                        if(subNodeHP[node]==0){
                            hp1Count++;
                        }
                        else{
                            hp2Count++;
                        }
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
            bestEdgeConnect[refEnd] = refStart;
            bestEdgeConnect[altEnd] = altStart;
            connect = true;
        }
        else if ( selectConnect == 2 ){
            PosAllele refStart = std::make_pair((*breakNode).first, 2);
            PosAllele altStart = std::make_pair((*breakNode).first, 1);  
            PosAllele refEnd   = std::make_pair((*nodeIter).first,  1);
            PosAllele altEnd   = std::make_pair((*nodeIter).first,  2);
            bestEdgeConnect[refEnd] = refStart;
            bestEdgeConnect[altEnd] = altStart;
            connect = true;
        }
    }
    return connect;
}


void VairiantGraph::storeResultPath(){
    // Initialize the parameters of the stored results
    posAppear.clear();
    blockStart.clear();
    bkResult.clear();
    subNodeHP.clear();
    blockVec.clear();

    // loop all best edge, try to restore the full connect path
    for(std::map<PosAllele,PosAllele>::iterator edgeIter = bestEdgeConnect.begin() ; edgeIter != bestEdgeConnect.end() ; edgeIter++ ){
        PosAllele A = (*edgeIter).first;
        PosAllele B = (*edgeIter).second;
        
        PosAllele aOtherSide = std::make_pair(A.first, (A.second != 1 ? 1 : 2));
        PosAllele bOtherSide = std::make_pair(B.first, (B.second != 1 ? 1 : 2));

        // in order to connect block path
        std::map<PosAllele,int>::iterator aBlockIter = bkResult.find(A);
        std::map<PosAllele,int>::iterator bBlockIter = bkResult.find(B);

        std::map<int,int>::iterator aAppearIter = posAppear.find(A.first);
        std::map<int,int>::iterator bAppearIter = posAppear.find(B.first);

        std::map<int,int>::iterator aStartIter = blockStart.find(A.first);
        std::map<int,int>::iterator bStartIter = blockStart.find(B.first);
        
        int minBlockStart = std::min(A.first+1, B.first+1);

        if( aStartIter != blockStart.end() )
            minBlockStart = std::min(minBlockStart, (*aStartIter).second );
        if( bStartIter != blockStart.end() )
            minBlockStart = std::min(minBlockStart, (*bStartIter).second );

        // two node doesn't appear, block start
        if( aBlockIter == bkResult.end() && bBlockIter == bkResult.end() ){
            bkResult[A] = minBlockStart;
            bkResult[B] = minBlockStart;
            // add HP tag
            if( aAppearIter == posAppear.end() && bAppearIter == posAppear.end() ){
                subNodeHP[A] = 0;
                subNodeHP[B] = 0;
                subNodeHP[aOtherSide] = 1;
                subNodeHP[bOtherSide] = 1;
            }
        }
        // node A doesn't appear
        else if( aBlockIter == bkResult.end() && bBlockIter != bkResult.end() ){
            bkResult[A] = std::min(minBlockStart, bkResult[B]);
            bkResult[aOtherSide] = std::min(minBlockStart, bkResult[B]);
            
            subNodeHP[A] = subNodeHP[B];
            subNodeHP[aOtherSide] = ( subNodeHP[A] != 0 ? 0 : 1);//subNodeHP[bOtherSide];
        }
        // node B doesn't appear
        else if( aBlockIter != bkResult.end() && bBlockIter == bkResult.end() ){
            bkResult[B] = std::min(minBlockStart, bkResult[A]);
            bkResult[bOtherSide] = std::min(minBlockStart, bkResult[A]);
            
            subNodeHP[B] = subNodeHP[A];
            subNodeHP[bOtherSide] = ( subNodeHP[B] != 0 ? 0 : 1);//subNodeHP[aOtherSide];
        }
        else{
            // For diploid,
            // In most cases, there will be two sub-edges between every two positions.
            // And the phasing of two positions has been completed on the first sub-edge.
        }

        posAppear[A.first] = 1;
        posAppear[B.first] = 1;
        
        blockStart[A.first] = minBlockStart;
        blockStart[B.first] = minBlockStart;
        
        if( aAppearIter == posAppear.end())
            blockVec[minBlockStart].push_back(A.first);
        if( bAppearIter == posAppear.end())
            blockVec[minBlockStart].push_back(B.first);
    }
}


VairiantGraph::VairiantGraph(std::string &in_ref, PhasingParameters &in_params){
    params=&in_params;
    ref=&in_ref;
}

VairiantGraph::~VairiantGraph(){
}

void VairiantGraph::addEdge(std::vector<ReadVariant> &in_readVariant){
    readVariant = &in_readVariant;
    // iter all read
    for(std::vector<ReadVariant>::iterator readIter = in_readVariant.begin() ; readIter != in_readVariant.end() ; readIter++ ){
        // assign SV quality
        for( auto variant : (*readIter).variantVec ){
            if( variant.quality == -1 ){
                svPosition[variant.position] = 1;
                if( variant.allele == 1 ){
                    // SV calling
                    variant.quality = 60; 
                }
                else{
                    // assume ref
                    variant.quality = 30;
                }
            }
            nodeInfo[variant.position][(*readIter).read_name] = variant.quality;
        }
        // iter all pair of snp and construct initial graph
        std::vector<Variant>::iterator variant1Iter = (*readIter).variantVec.begin();
        std::vector<Variant>::iterator variant2Iter = std::next(variant1Iter,1);
        while(variant1Iter != (*readIter).variantVec.end() && variant2Iter != (*readIter).variantVec.end() ){
            // create new edge if not exist
            std::map<int,VariantEdge*>::iterator posIter = edgeList.find((*variant1Iter).position);
            if( posIter == edgeList.end() )
                edgeList[(*variant1Iter).position] = new VariantEdge((*variant1Iter).position);

            // add edge process
            for(int nextNode = 0 ; nextNode < params->connectAdjacent; nextNode++){
                // this allele support ref
                if( (*variant1Iter).allele == 0 )
                    edgeList[(*variant1Iter).position]->ref->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).read_name);
                // this allele support alt
                if( (*variant1Iter).allele == 1 )
                    edgeList[(*variant1Iter).position]->alt->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).read_name);
                
                // next snp
                variant2Iter++;
                if( variant2Iter == (*readIter).variantVec.end() ){
                    break;
                }
            }

            variant1Iter++;
            variant2Iter = std::next(variant1Iter,1);
        }
    }
} 

void VairiantGraph::readCorrection(){
    //std::map<int,int> inconsistentSnpMap;
    std::map<int,int> snpAlleleCountMap;
    //std::map<int,int> deleteSnpMap;
    std::map<std::string,int> readHP;

    // iter all read, determine the haplotype of the read and mark abnormal SNPs
    for(std::vector<ReadVariant>::iterator readIter = (*readVariant).begin() ; readIter != (*readVariant).end() ; readIter++ ){
        double refCount = 0;
        double altCount = 0;
        // loop all variant 
        for( auto variant : (*readIter).variantVec ){
            PosAllele node = std::make_pair( variant.position , variant.allele+1);
            std::map<PosAllele,int>::iterator nodePS = bkResult.find(node);

            if( nodePS != bkResult.end() ){
                snpAlleleCountMap[variant.position]++;
                
                if( bkResult[node] != 0 ){
                    if(subNodeHP[node]==0)refCount++;
                    else altCount++;
                }
            }
        }

        // mark abnormal SNPs. the result of tagging is determined by phased alleles
        if( std::max(refCount,altCount)/(refCount+altCount) >= 0.6 ){
            int belongHP = ( refCount > altCount ? 0 : 1 );
            readHP[(*readIter).read_name] = belongHP;
        }
        else{
            readHP[(*readIter).read_name] = -1;
        }
    }

    std::map<int,std::map<int,std::map<int,int>>> hpAlleleCountMap;

    // count each haplotype's allele on SNP locus
    for(auto readIter = (*readVariant).begin() ; readIter != (*readVariant).end() ; readIter++ ){
        
        auto hpIter = readHP.find((*readIter).read_name);
        if( hpIter == readHP.end() )
            continue;
        if( readHP[(*readIter).read_name] == -1 )
            continue;
        
        int rHP = readHP[(*readIter).read_name];
        for(auto variantIter = (*readIter).variantVec.begin() ; variantIter != (*readIter).variantVec.end() ; variantIter++ ){
            if( (*variantIter).allele == 0 || (*variantIter).allele == 1){
                int allele = (*variantIter).allele;
                hpAlleleCountMap[rHP][(*variantIter).position][allele]++;
            }
        }
    }
  
    subNodeHP.clear();
    
    std::map<int,std::map<int,int>> hpAllele;
    // reassign allele result
    for(auto nodeIter = nodeInfo.begin() ; nodeIter != nodeInfo.end() ; nodeIter++ ){
        int position = (*nodeIter).first;
        PosAllele A = std::make_pair(position, 1);
        PosAllele aOtherSide = std::make_pair(position, 2);
        
        double hp1R = hpAlleleCountMap[0][(*nodeIter).first][0];
        double hp1A = hpAlleleCountMap[0][(*nodeIter).first][1];
        double hp2R = hpAlleleCountMap[1][(*nodeIter).first][0];
        double hp2A = hpAlleleCountMap[1][(*nodeIter).first][1];
        double hp1ConsistentRatio = std::max(hp1R,hp1A)/(hp1R + hp1A);
        double hp2ConsistentRatio = std::max(hp2R,hp2A)/(hp2R + hp2A);
        double maxAlleleRatio = std::max(hp1R + hp2R, hp1A + hp2A)/(hp1R + hp2R + hp1A + hp2A);

        if(hp1ConsistentRatio >= params->alleleConsistentRatio){
            // decide hp1 allele. at least two identical alleles are required to determine haplotype
            if( hp1R > hp1A && hp1R >= 2 )
                hpAllele[1][position] = 0;
            else if( hp1R < hp1A && hp1A >= 2 )
                hpAllele[1][position] = 1;
            else
                hpAllele[1][position] = -1;
        }
        else{
            hpAllele[1][position] = -1;
        }
        
        if(hp2ConsistentRatio >= params->alleleConsistentRatio){
            // decide hp2 allele. at least two identical alleles are required to determine haplotype
            if( hp2R > hp2A && hp2R >= 2 )
                hpAllele[2][position] = 0;
            else if( hp2R < hp2A && hp2A >= 2 )
                hpAllele[2][position] = 1;
            else
                hpAllele[2][position] = -1;
        }
        else{
            hpAllele[2][position] = -1; 
        }

        if( maxAlleleRatio > params->maxAlleleRatio ){
            if( ( hp1R > hp1A && hp2R > hp2A ) || ( hp1R < hp1A && hp2R < hp2A ) ){
                hpAllele[1][position] = -1; 
                hpAllele[2][position] = -1; 
            }
        }
        
        if( hpAllele[1][position] == hpAllele[2][position] ){
            hpAllele[1][position] = -1;
            hpAllele[2][position] = -1;
        }
        else if( hpAllele[1][position] == -1 && hpAllele[2][position] != -1){
            hpAllele[1][position] = (hpAllele[2][position] == 0 ? 1 : 0);
        }
        else if( hpAllele[2][position] == -1 && hpAllele[1][position] != -1){
            hpAllele[2][position] = (hpAllele[1][position] == 0 ? 1 : 0);
        }
        
        if( hpAllele[1][position] == 0 || hpAllele[2][position] == 1 ){
            subNodeHP[A] = 0;
            subNodeHP[aOtherSide] = 1;
        }
        else if( hpAllele[1][position] == 1 || hpAllele[2][position] == 0 ){
            subNodeHP[A] = 1;
            subNodeHP[aOtherSide] = 0;
        }
        else{
            bkResult.erase(A);
            bkResult.erase(aOtherSide);
        }  
    }
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
    for( std::map<int,ReadBase>::iterator nodeIter = nodeInfo.begin() ; nodeIter != nodeInfo.end() ; nodeIter++ ){
        
        PhasingElement tmp;
        
        PosAllele ref = std::make_pair( (*nodeIter).first , 1);
        PosAllele alt = std::make_pair( (*nodeIter).first , 2);
        
        std::map<PosAllele,int>::iterator psRefIter = bkResult.find(ref);
        std::map<PosAllele,int>::iterator psAltIter = bkResult.find(alt);
        
        if( psRefIter != bkResult.end() || psAltIter != bkResult.end() ){
            if( psRefIter != bkResult.end() )
                tmp.block = (*psRefIter).second;
            else
                tmp.block = (*psAltIter).second;
            tmp.RAstatus = std::to_string(subNodeHP[ref]) + "|" + std::to_string(subNodeHP[alt]);
        }
        else
            continue;
        
        if( tmp.block != 0){
            std::string key = chrName + "_" + std::to_string( (*nodeIter).first );
            result[key] = tmp;
        }
    }
}

int VairiantGraph::totalNode(){
    return nodeInfo.size();
}

void VairiantGraph::phasingProcess(){

    this->edgeConnectResult();
    this->storeResultPath();
    /*
    while(true){
        //bool connect = this->connectBlockByCommonRead();
        bool connect = this->blockPhaseing();
        this->storeResultPath();

        if(!connect){
            break;
        }
    }
    */
    this->readCorrection();  
}


  