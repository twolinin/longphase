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

std::vector<std::string> SubEdge::getSupportRead(int breakNodePosition){

    std::vector<std::string> readVec;
    for(std::map<int, std::vector<std::string> >::iterator edgeIter = refRead.begin() ; edgeIter != refRead.end() ; edgeIter++ ){
        for(std::vector<std::string>::iterator readIter = (*edgeIter).second.begin() ; readIter != (*edgeIter).second.end() ; readIter++ ){
            readVec.push_back((*readIter));
        }
    }
    
    for(std::map<int, std::vector<std::string> >::iterator edgeIter = altRead.begin() ; edgeIter != altRead.end() ; edgeIter++ ){
        for(std::vector<std::string>::iterator readIter = (*edgeIter).second.begin() ; readIter != (*edgeIter).second.end() ; readIter++ ){
            readVec.push_back((*readIter));
        }
    }
    return readVec;
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
std::pair<PosAllele,PosAllele> VariantEdge::findBestEdgePair(int targetPos, bool isONT, double diffRatioThreshold, double svDiffRatioThreshold, bool containSV){
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
        if(rr!=0 && aa !=0){
            rr+=1;
            aa+=1;
        }
        if(ra!=0 && ar !=0){
            ra+=1;
            ar+=1;
        }
    }

    double readsThreshold = (double)std::min((rr+aa),(ar+ra)) / (double)std::max((rr+aa),(ar+ra));

    if( readsThreshold > diffRatioThreshold && !containSV){
        // SNP-SNP
        // The number of RR and RA is similar, and it is easily affected by sequence errors.
    }
    else if( readsThreshold > svDiffRatioThreshold && containSV){
        // SNP-SV or SV-SNP
        // The number of RR and RA is similar, and it is easily affected by sequence errors.
    }
    else if( rr + aa > ra + ar ){
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

//VairiantGrpah
void VairiantGrpah::initialResult(){
    
    // current snp, haplotype (1 or 2), support snp
    std::map<int, std::map<int,std::vector<int> > > hpCountMap;
    // current snp, result haplotype (1 or 2)
    std::map<int,int> hpResult;
    // current snp, result haplotype (1 or 2)
    std::map<int,int> inconsistentSnpMap;
    // < block start, <snp in this block> >
    std::map<int,std::vector<int> > phasedBlocks;
    int blockStart;
    
    // Visit all position and assign SNPs to haplotype
    // Only one of the two alleles needs to be used for each SNP
    for(std::map<int,ReadBase>::iterator nodeIter = nodeInfo.begin() ; nodeIter != nodeInfo.end() ; nodeIter++ ){
        // check next position
        std::map<int,ReadBase>::iterator nextNodeIter = std::next(nodeIter, 1);
        if( nextNodeIter == nodeInfo.end() )
             break;
        
        int currPos = (*nodeIter).first;
        int nextPos = (*nextNodeIter).first;

        // check distance between two position 
        if(std::abs(nextPos-currPos) > params->distance ){
            continue;
        }
        
        // get the number of HP1 and HP2 supported reference allele
        int h1 = hpCountMap[currPos][1].size();
        int h2 = hpCountMap[currPos][2].size();
        
        // check and record inconsistent snp
        if( (h1==1 || h2==1) && h1 + h2 > 3 ){
            // if number of HP1 supported is 1, this supported snp is inconsistent
            int inconsistentSnp = hpCountMap[currPos][ (h1==1 ? 1 : 2) ][0];
            inconsistentSnpMap[inconsistentSnp]++;
        }
        
        // new block, set this position as block start 
        if( h1 == 0 && h2 == 0 ){
            if( hpCountMap[nextPos][1].size() !=0 || hpCountMap[nextPos][2].size() !=0 ){
                continue;
            }
            blockStart=currPos;
            phasedBlocks[blockStart].push_back(currPos);
            hpResult[currPos] = 1;
        }
        // ignore poorly classified SNPs
        else if( (h1==2&&h2==3) || (h1==3&&h2==2) || (h1==2&&h2==1) || (h1==1&&h2==2) ){
            hpResult[currPos] = 0;
            continue;
        }
        else{
            if( h1 > h2 || h1 < h2 ){
                int currHP = ( h1 > h2 ? 1 : 2 );
                hpResult[currPos] = currHP;
                phasedBlocks[blockStart].push_back(currPos);
            }
            else{
                // equal weight node, go back to judge result
                std::map<int,ReadBase>::iterator preStartIter = std::prev(nodeIter, 1);
                int preStart = (*preStartIter).first;
                int preNext  = (*nodeIter).first;
                // check previos SNP in edge list
                std::map<int,VariantEdge*>::iterator edgeIter = edgeList.find( preStart );
                if( edgeIter==edgeList.end() ){
                    hpResult[currPos] = 0;
                    continue;
                }
                // just consider reads from the previous SNP and currnt SNP 
                std::pair<PosAllele,PosAllele> preResult = (*edgeIter).second->findBestEdgePair(preNext, params->isONT, 1, 1, false);

                // greedy RR or RA by number of read 
                if( preResult.first.second == 1 ){
                    hpResult[currPos] = hpResult[preStart];
                }
                else if( preResult.first.second == 2 ){
                    hpResult[currPos] = (hpResult[preStart] == 1 ? 2 : 1);
                }
                else{
                    // RR and RA are equal read, judge by base quality
                    PosAllele refEnd = std::make_pair(preNext, 1);
                    PosAllele altEnd = std::make_pair(preNext, 2);
                    
                    // get four edge quality
                    int rrAvgQuality = edgeList[preStart]->ref->getAvgQuality(refEnd);
                    int raAvgQuality = edgeList[preStart]->ref->getAvgQuality(altEnd);
                    int arAvgQuality = edgeList[preStart]->alt->getAvgQuality(refEnd);
                    int aaAvgQuality = edgeList[preStart]->alt->getAvgQuality(altEnd);
                    
                    if( rrAvgQuality + aaAvgQuality > raAvgQuality + arAvgQuality )
                        hpResult[currPos] = hpResult[preStart];
                    else
                        hpResult[currPos] = (hpResult[preStart] == 1 ? 2 : 1);
                }
            }
        }
        
        // Check if there is no edge from current node
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList.find( currPos );
        if( edgeIter==edgeList.end() )
            continue;
        
        // check connect between surrent SNP and next 5 SNPs
        for(int i = 0 ; i < 5 ; i++ ){
            // consider reads from the currnt SNP and the next (i+1)'s SNP
            std::pair<PosAllele,PosAllele> tmp = (*edgeIter).second->findBestEdgePair((*nextNodeIter).first, params->isONT, 1, 1, false);
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
            }
            nextNodeIter++;
            if( nextNodeIter == nodeInfo.end() )
                break;
        }
    }

    // Correction for inconsistent SNPs
    std::vector<std::string>::iterator edgeIter = dotResult.begin();
    for(std::map<int,int>::iterator snpIter = inconsistentSnpMap.begin();snpIter != inconsistentSnpMap.end();snpIter++){
        if((*snpIter).second>=4){
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
        if( (*blockIter).second.size()<=1 )
            continue;
        
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

std::map<std::string,int> VairiantGrpah::getBlockRead(std::pair<int,std::vector<int> > currentBlockVec, std::map<std::string,int> &readQuality, BlockRead &totalRead , int sampleNum){
    //BlockRead hp1Read;
    //BlockRead hp2Read;
    BlockRead currentTotalRead;
    std::map<std::string,int> readTag;
    std::map<std::string,int> readTmpQuality;
    
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
        
        // get ref reads
        std::vector<std::string> rReadVec = (*edgeIter).second->ref->getSupportRead((*posIter));
        // tag REF haplotype read
        for(std::vector<std::string>::iterator readIter = rReadVec.begin() ; readIter != rReadVec.end() ; readIter++ ){
            // Cumulative read quality
            std::map<std::string,int>::iterator readQuaIter = readTmpQuality.find((*readIter));
            if( readQuaIter == readTmpQuality.end() ){
                readTmpQuality[(*readIter)] = nodeInfo[(*posIter)][(*readIter)];
            }
            else{
                readTmpQuality[(*readIter)] += nodeInfo[(*posIter)][(*readIter)];
            }
            currentTotalRead.recordRead((*readIter));
            totalRead.recordRead((*readIter));
        }
        
        // get alt reads
        std::vector<std::string> aReadVec = (*edgeIter).second->alt->getSupportRead((*posIter));
        // tag ALT haplotype read    
        for(std::vector<std::string>::iterator readIter = aReadVec.begin() ; readIter != aReadVec.end() ; readIter++ ){ 
            // Cumulative read quality
            std::map<std::string,int>::iterator readQuaIter = readTmpQuality.find((*readIter));
            if( readQuaIter == readTmpQuality.end() ){
                readTmpQuality[(*readIter)] = nodeInfo[(*posIter)][(*readIter)] * (-1);
            }
            else{
                readTmpQuality[(*readIter)] -= nodeInfo[(*posIter)][(*readIter)];
            }       
            currentTotalRead.recordRead((*readIter));
            totalRead.recordRead((*readIter));
        }
    }   
    // Divide read into two different haplotypes
    for( std::map<std::string,int>::iterator readIter = currentTotalRead.readVec.begin() ; readIter != currentTotalRead.readVec.end() ; readIter++ ){
        if( readTmpQuality[(*readIter).first] > 0 ){
            readTag[(*readIter).first] = 0;
        }
        else if( readTmpQuality[(*readIter).first] < 0){
            readTag[(*readIter).first] = 1;
        }
        else{
        }
        readQuality[(*readIter).first] = std::abs(readTmpQuality[(*readIter).first]);
    }
    return readTag;
}

bool VairiantGrpah::connectBlockByCommonRead(int nextBlcok, double diffRatioThreshold){
    bool connect = false;
    // loop all block, prepare to create new edges between two block
    // PS number, < position vector >
    for(std::map<int,std::vector<int> >::iterator blockVecIter = blockVec.begin() ; blockVecIter != blockVec.end() ; blockVecIter++ ){
        // Number of SNPs checked
        int sampleNum = 15;
        // get next block
        std::map<int,std::vector<int> >::iterator nextBlockIter = std::next(blockVecIter, nextBlcok);
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
        std::map<int,std::vector<int> >::iterator midBlockIter = std::next(blockVecIter, 1);

        if(nextBlcok==2){
            // Some small blocks will be in the middle of two large blocks.
            // This small block looks like an island.
            // if the middle block is very large, it should not be crossed over
            std::map<int,std::vector<int> >::iterator midBlockIter = std::next(blockVecIter, 1);
            int midBlockStart = (*midBlockIter).second[0];
            int midBlockEnd   = (*midBlockIter).second[(*midBlockIter).second.size()-1];
            // middle block size threshold
            if( midBlockEnd - midBlockStart >= params->islandBlockLength ){
                continue;
            }
        }

        BlockRead totalRead;
        std::map<std::string,int> frontReadQuality;
        std::map<std::string,int> backReadQuality;
        std::map<std::string,int> frontReadtag = this->getBlockRead((*blockVecIter), frontReadQuality, totalRead, sampleNum);
        std::map<std::string,int> backReadtag  = this->getBlockRead((*nextBlockIter), backReadQuality, totalRead, sampleNum);
        int rr=0;
        int ra=0;
        int ar=0;
        int aa=0;
        int rrTmpQyality=0;
        int raTmpQyality=0;
        int arTmpQyality=0;
        int aaTmpQyality=0;
        int selectConnect = -1;
        
        // loop all block read, judge each read connect result
        for( std::map<std::string,int>::iterator readIter = totalRead.readVec.begin() ; readIter != totalRead.readVec.end() ; readIter++ ){
            std::map<std::string,int>::iterator frontReadIter = frontReadtag.find((*readIter).first);
            std::map<std::string,int>::iterator backReadIter  = backReadtag.find((*readIter).first);
            // this read exists in two blocks
            if( frontReadIter != frontReadtag.end() && backReadIter != backReadtag.end() ){
                if( (*frontReadIter).second == 0 && (*backReadIter).second == 0 ){
                    rrTmpQyality += frontReadQuality[(*readIter).first] + backReadQuality[(*readIter).first];
                    rr++;
                }
                else if( (*frontReadIter).second == 0 && (*backReadIter).second == 1 ){
                    raTmpQyality += frontReadQuality[(*readIter).first] + backReadQuality[(*readIter).first];
                    ra++;
                }
                else if( (*frontReadIter).second == 1 && (*backReadIter).second == 0 ){
                    arTmpQyality += frontReadQuality[(*readIter).first] + backReadQuality[(*readIter).first];
                    ar++;
                } 
                else if( (*frontReadIter).second == 1 && (*backReadIter).second == 1 ){
                    aaTmpQyality += frontReadQuality[(*readIter).first] + backReadQuality[(*readIter).first];
                    aa++;
                }
            }
        }
        // two sides in a combination with higher credibility
        if( rr + aa == ra + ar ){
            if(rr!=0 && aa !=0){
                rr+=1;
                aa+=1;
            }
            if(ra!=0 && ar !=0){
                ra+=1;
                ar+=1;
            }
        }

        if( rr + aa > ra + ar ){
            selectConnect = 1;
        }
        else if(rr + aa < ra + ar){
            selectConnect = 2;
        }
        else{
            // weight(aa + rr) == weight(ar+ ra)
            // the edges are equal but the quality may be different
            if( rr + aa + ra + ar !=0 ){
                if( rrTmpQyality + aaTmpQyality > raTmpQyality + arTmpQyality ){
                    selectConnect = 1;
                }
                else if(rrTmpQyality + aaTmpQyality < raTmpQyality + arTmpQyality){
                    selectConnect = 2;
                }
                else{
                }
            }
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

void VairiantGrpah::findResultPath(){
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


VairiantGrpah::VairiantGrpah(std::string &in_ref, PhasingParameters &in_params){
    params=&in_params;
    ref=&in_ref;
}

VairiantGrpah::~VairiantGrpah(){
}

void VairiantGrpah::addEdge(std::vector<ReadVariant> &in_readVariant){
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
            for(int nextNode = 0 ; nextNode < 5; nextNode++){
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

void VairiantGrpah::readTest(){
    std::map<int,int> inconsistentSnpMap;
    std::map<int,int> snpAlleleCountMap;
    std::map<int,int> deleteSnpMap;

    // iter all read, determine which haplotype each read belongs to and mark abnormal SNPs
    for(std::vector<ReadVariant>::iterator readIter = (*readVariant).begin() ; readIter != (*readVariant).end() ; readIter++ ){
        double refCount = 0;
        double altCount = 0;
        // loop all variant 
        for( auto variant : (*readIter).variantVec ){
            PosAllele node = std::make_pair( variant.position , variant.allele+1);
            std::map<PosAllele,int>::iterator nodePS = bkResult.find(node);

            if( nodePS != bkResult.end() ){
                snpAlleleCountMap[variant.position]++;
                
                if(subNodeHP[node]==0)refCount++;
                else altCount++;
            }
        }

        // mark abnormal SNPs
        if( std::max(refCount,altCount)/(refCount+altCount) >= 0.6 ){
            int belongHP = ( refCount > altCount ? 0 : 1 );
            // loop all variant
            for( auto variant : (*readIter).variantVec ){
                PosAllele node = std::make_pair( variant.position , variant.allele+1);
                std::map<PosAllele,int>::iterator nodePS = bkResult.find(node);
                // count inconsistent snp
                if( nodePS != bkResult.end() ){
                    if( subNodeHP[node] != belongHP ){
                        inconsistentSnpMap[variant.position]++;
                    }
                }
            }
        }
    }
    // identify abnormal SNPs
    for(auto snpIter = inconsistentSnpMap.begin() ; snpIter != inconsistentSnpMap.end() ; snpIter++ ){
        double totalCount = snpAlleleCountMap[(*snpIter).first];
        double inConsisentCount = (*snpIter).second;
        double consisentCount   = totalCount - inConsisentCount;
        double inConsisentRatio = std::min(inConsisentCount,consisentCount)/totalCount;

        if( inConsisentRatio >= 0.4 ){
            deleteSnpMap[(*snpIter).first]=1;
        }
    }

    // erase abnormal SNPs
    for( auto snpIter = deleteSnpMap.begin() ; snpIter != deleteSnpMap.end() ; snpIter++ ){
        auto nodeIter = nodeInfo.find((*snpIter).first);
        if( nodeIter != nodeInfo.end() )
            nodeInfo.erase(nodeIter);
        
        auto svIter = svPosition.find((*snpIter).first);
        if( svIter != svPosition.end() )
            svPosition.erase(svIter);
    }

    // filter abnormal SNPs
    for(auto readIter = (*readVariant).begin() ; readIter != (*readVariant).end() ; readIter++ ){
        auto variantIter = (*readIter).variantVec.begin();
        while( variantIter != (*readIter).variantVec.end() ){
            auto delIter = deleteSnpMap.find((*variantIter).position);
            if( delIter != deleteSnpMap.end() ){
                if( variantIter == (*readIter).variantVec.begin() ){
                    (*readIter).variantVec.erase(variantIter);
                    continue;
                }
                else{
                    (*readIter).variantVec.erase(variantIter);
                    variantIter--;
                }
            }
            variantIter++;
        }
    }

    edgeList.clear();
    bestEdgeConnect.clear();
    dotResult.clear();
    
    // reconnect construct edge
    // iter all pair of snp and construct initial graph
    for(auto readIter = (*readVariant).begin() ; readIter != (*readVariant).end() ; readIter++ ){
        std::vector<Variant>::iterator variant1Iter = (*readIter).variantVec.begin();
        std::vector<Variant>::iterator variant2Iter = std::next(variant1Iter,1);
        while(variant1Iter != (*readIter).variantVec.end() && variant2Iter != (*readIter).variantVec.end() ){
            // create new edge if not exist
            std::map<int,VariantEdge*>::iterator posIter = edgeList.find((*variant1Iter).position);
            if( posIter == edgeList.end() )
                edgeList[(*variant1Iter).position] = new VariantEdge((*variant1Iter).position);

            // add edge process
            for(int nextNode = 0 ; nextNode < 5; nextNode++){
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

void VairiantGrpah::writingDotFile(std::string dotPrefix){
    
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

void VairiantGrpah::exportResult(std::string chrName, PhasingResult &result){
    
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
        
        std::string key = chrName + "_" + std::to_string( (*nodeIter).first );
        result[key] = tmp;
    }
}

int VairiantGrpah::totalNode(){
    return nodeInfo.size();
}

void VairiantGrpah::phasingProcess(){

    for( int round = 1 ; round <= 2 ; round++ ){
        this->initialResult();
        this->findResultPath();
      
        for(int i = 1 ; i <= params->crossBlock ; i++ ){
            while(true){
                bool connect = this->connectBlockByCommonRead(i, params->blockReadThreshold);
                this->findResultPath();

                if(!connect){
                    break;
                }
            }
        }
        
        if( round == 1 )
            this->readTest();
    }
}


  