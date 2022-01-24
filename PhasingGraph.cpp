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
            return -1;
        else
            return refQuality[targetPos.first];
    }
    // target is Alt allele
    if( targetPos.second == 2 ){
        std::map<int, int>::iterator qIter = altQuality.find(targetPos.first);
        if( qIter == altQuality.end() )
            return -1;
        else
            return altQuality[targetPos.first];
    }
    return -1;
}

int SubEdge::getAvgQuality(PosAllele targetPos){
    // target is Ref allele
    if( targetPos.second == 1 ){
        std::map<int, int>::iterator qIter = refQuality.find(targetPos.first);
        if( qIter == refQuality.end() )
            return -1;
        else
            return refQuality[targetPos.first]/refReadCount[targetPos.first];
    }
    // target is Alt allele
    if( targetPos.second == 2 ){
        std::map<int, int>::iterator qIter = altQuality.find(targetPos.first);
        if( qIter == altQuality.end() )
            return -1;
        else
            return altQuality[targetPos.first]/altReadCount[targetPos.first];
    }
    return -1;
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



void VairiantGrpah::findBestEdgeList(){
    // Check whether there are edges in adjacent position
    // and select the edge with the largest weight combination 
    // as the best choice for the subgraph

    // visit all nodes
    for(std::map<int,ReadBase>::iterator nodeIter = nodeInfo.begin() ; nodeIter != nodeInfo.end() ; nodeIter++ ){
        // get next node
        std::map<int,ReadBase>::iterator nextNodeIter = std::next(nodeIter, 1);
        if( nextNodeIter == nodeInfo.end() )
             break;
        
        int currPos = (*nodeIter).first;
        int nextPos = (*nextNodeIter).first;

        // check distance between two position 
        if(std::abs(nextPos-currPos) > params->distance ){
            continue;
        }
        
        // Check if there is an edge from current node
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList.find( currPos );
        if( edgeIter==edgeList.end() ){
            continue;
        }
        
        std::map<int,int>::iterator currSVIter = svPosition.find(currPos);
        std::map<int,int>::iterator nextSVIter = svPosition.find(nextPos);
        
        bool currSV = ( currSVIter != svPosition.end() );
        bool nextSV = ( nextSVIter != svPosition.end() );
        bool containSV = (currSV || nextSV);
        
        // find best edge pair from current position to next position
        std::pair<PosAllele,PosAllele> EdgePair = (*edgeIter).second->findBestEdgePair(nextPos, params->isONT, params->readsThreshold, params->svReadsThreshold, containSV);

        // restore the position and allele of edge pair
        PosAllele refStart = std::make_pair(currPos, 1);
        PosAllele refEnd   = EdgePair.first;
        PosAllele altStart = std::make_pair(currPos, 2);
        PosAllele altEnd   = EdgePair.second;

        // -1 : ref or allele at this position didn't find the best edge to the next
        if( refEnd.second != -1 )
            bestEdgeConnect[refEnd] = refStart;
        if( altEnd.second != -1 )
            bestEdgeConnect[altEnd] = altStart;
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

bool VairiantGrpah::connectBlockByReadName(int nextBlcok, double diffRatioThreshold, bool blockPhasing, bool final){
    bool connect = false;
    std::map<int,std::vector<int> >::reverse_iterator last = blockVec.rbegin() ;
    
    // loop all block, prepare to create new edges between two block
    // PS number, < position vector >
    for(std::map<int,std::vector<int> >::iterator blockVecIter = blockVec.begin() ; blockVecIter != blockVec.end() ; blockVecIter++ ){
        // Number of SNPs checked
        int sampleNum = 1 ;   
        if(blockPhasing)
            sampleNum = 15;
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
        if( edgeIter==edgeList.end() && blockPhasing){
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
            if( midBlockEnd - midBlockStart >= params->islandBlockLength && blockPhasing){
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
        if( blockPhasing && rr + aa == ra + ar ){
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
            if( blockPhasing && rr + aa + ra + ar !=0 ){
                if( rrTmpQyality + aaTmpQyality > raTmpQyality + arTmpQyality ){
                    selectConnect = 1;
                }
                else if(rrTmpQyality + aaTmpQyality < raTmpQyality + arTmpQyality){
                    selectConnect = 2;
                }
                else{
                }
            }
            // The number of reads is the same in the RR and RA of the two blocks.
            // Finally, return to the read between the two SNPs to make a decision.
            else if(final){
                // get breakpoint SNP
                std::map<int,VariantEdge*>::iterator edgeIter = edgeList.find( (*breakNode).first );
                // get next SNP
                std::map<int,ReadBase>::iterator currNode = breakNode;
                std::map<int,ReadBase>::iterator nextNode = std::next(breakNode, 1);

                while( currNode != nodeInfo.end() && nextNode != nodeInfo.end() && (*currNode).first < (*nodeIter).first ){
                    std::map<int,VariantEdge*>::iterator edgeIter = edgeList.find( (*currNode).first );
                    
                    if(edgeIter != edgeList.end() ){
                        std::pair<PosAllele,PosAllele> EdgePair = (*edgeIter).second->findBestEdgePair((*nextNode).first, params->isONT, diffRatioThreshold, 0, false);
                        // restore the position and allele of edge pair
                        PosAllele refStart = std::make_pair((*currNode).first, 1);
                        PosAllele refEnd   = EdgePair.first;
                        PosAllele altStart = std::make_pair((*currNode).first, 2);
                        PosAllele altEnd   = EdgePair.second;
                        
                        // -1 : ref or allele at this position didn't find the best edge to the next
                        if( refEnd.second != -1 ){
                            bestEdgeConnect[refEnd] = refStart;
                        }
                        if( altEnd.second != -1 ){
                            bestEdgeConnect[altEnd] = altStart;
                        }
                    }
                    currNode++;
                    nextNode++;
                }
                continue;
            }
        }

        double readsThreshold = (double)std::min((rr+aa),(ar+ra)) / (double)std::max((rr+aa),(ar+ra));

        // Assign read to REF haplotype or ALT haplotype
        if( readsThreshold > diffRatioThreshold && !blockPhasing ){
            std::map<int,ReadBase>::iterator nextNode = std::next(breakNode, 1);

            if( rr+aa+ra+ar > 0 && !blockPhasing && (*nextNode).first == (*nodeIter).first ){
                
                PosAllele refEnd = std::make_pair((*nodeIter).first, 1);
                PosAllele altEnd = std::make_pair((*nodeIter).first, 2);
                
                // get four edge quality
                int rrAvgQuality = edgeList[(*breakNode).first]->ref->getAvgQuality(refEnd);
                int raAvgQuality = edgeList[(*breakNode).first]->ref->getAvgQuality(altEnd);
                int arAvgQuality = edgeList[(*breakNode).first]->alt->getAvgQuality(refEnd);
                int aaAvgQuality = edgeList[(*breakNode).first]->alt->getAvgQuality(altEnd);
                double qualityThreshold = (double)std::min((rrAvgQuality+aaAvgQuality),(raAvgQuality+arAvgQuality)) / (double)std::max((rrAvgQuality+aaAvgQuality),(raAvgQuality+arAvgQuality));
                
                if( qualityThreshold > diffRatioThreshold ){
                }
                else if ( rrAvgQuality+aaAvgQuality > raAvgQuality+arAvgQuality ){
                    selectConnect = 1;
                }
                else if ( rrAvgQuality+aaAvgQuality < raAvgQuality+arAvgQuality ){
                    selectConnect = 2;
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

void VairiantGrpah::connectBlockByQuality(){
    // loop all block, prepare to create new edges between two block
    // PS number, < position vector >
    for(std::map<int,std::vector<int> >::iterator currentBlockIter = blockVec.begin() ; currentBlockIter != blockVec.end() ; currentBlockIter++ ){
        // The last SNP of this block
        std::vector<int>::reverse_iterator posIter = (*currentBlockIter).second.rbegin();
        // current SNP iter
        std::map<int,ReadBase>::iterator currentPosIter = nodeInfo.find((*posIter));

        // Check there have edges from current node to the next
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList.find( (*posIter) );
        if( edgeIter==edgeList.end() )
            continue;
        
        if( currentPosIter != nodeInfo.end() ){
            std::map<int,ReadBase>::iterator  nextPosIter = std::next(currentPosIter, 1);

            int currPos = (*posIter);
            int nextPos = (*nextPosIter).first;
              
            // check distance between two position 
            if(std::abs(nextPos-currPos) > params->distance )
                continue;
            
            PosAllele refEnd = std::make_pair((*nextPosIter).first, 1);
            PosAllele altEnd = std::make_pair((*nextPosIter).first, 2);
            // get four edge quality
            
            int rrQuality = edgeList[(*posIter)]->ref->getQuality(refEnd);
            int raQuality = edgeList[(*posIter)]->ref->getQuality(altEnd);
            int arQuality = edgeList[(*posIter)]->alt->getQuality(refEnd);
            int aaQuality = edgeList[(*posIter)]->alt->getQuality(altEnd);
                    
            if( rrQuality + aaQuality > raQuality + arQuality ){
                PosAllele refStart = std::make_pair((*posIter), 1);
                PosAllele altStart = std::make_pair((*posIter), 2);
                PosAllele refEnd   = std::make_pair((*nextPosIter).first,  1);
                PosAllele altEnd   = std::make_pair((*nextPosIter).first,  2);
                bestEdgeConnect[refEnd] = refStart;
                bestEdgeConnect[altEnd] = altStart;
            }
            else if( rrQuality + aaQuality < raQuality + arQuality){
                PosAllele refStart = std::make_pair((*posIter), 2);
                PosAllele altStart = std::make_pair((*posIter), 1);
                PosAllele refEnd   = std::make_pair((*nextPosIter).first,  1);
                PosAllele altEnd   = std::make_pair((*nextPosIter).first,  2);
                bestEdgeConnect[refEnd] = refStart;
                bestEdgeConnect[altEnd] = altStart;
            }

        }
    }
}

void VairiantGrpah::connectByQuality(){
    // visit all nodes
    for(std::map<int,ReadBase>::iterator nodeIter = nodeInfo.begin() ; nodeIter != nodeInfo.end() ; nodeIter++ ){
        // get next node
        std::map<int,ReadBase>::iterator nextNodeIter = std::next(nodeIter, 1);
        if( nextNodeIter == nodeInfo.end() )
             break;
        // Check if there is an edge from current node
        std::map<int,VariantEdge*>::iterator edgeIter = edgeList.find( (*nodeIter).first );
            if( edgeIter==edgeList.end() )
                continue;
                
        PosAllele refEnd = std::make_pair((*nextNodeIter).first, 1);
        PosAllele altEnd = std::make_pair((*nextNodeIter).first, 2);
        // get four edge quality
        int rrQuality = edgeList[(*nodeIter).first]->ref->getAvgQuality(refEnd);
        int raQuality = edgeList[(*nodeIter).first]->ref->getAvgQuality(altEnd);
        int arQuality = edgeList[(*nodeIter).first]->alt->getAvgQuality(refEnd);
        int aaQuality = edgeList[(*nodeIter).first]->alt->getAvgQuality(altEnd);

        if( rrQuality + aaQuality > raQuality + arQuality ){
            PosAllele refStart = std::make_pair((*nodeIter).first, 1);
            PosAllele altStart = std::make_pair((*nodeIter).first, 2);
            PosAllele refEnd   = std::make_pair((*nextNodeIter).first,  1);
            PosAllele altEnd   = std::make_pair((*nextNodeIter).first,  2);
            bestEdgeConnect[refEnd] = refStart;
            bestEdgeConnect[altEnd] = altStart;
        }
        else if( rrQuality + aaQuality < raQuality + arQuality){
            PosAllele refStart = std::make_pair((*nodeIter).first, 2);
            PosAllele altStart = std::make_pair((*nodeIter).first, 1);
            PosAllele refEnd   = std::make_pair((*nextNodeIter).first,  1);
            PosAllele altEnd   = std::make_pair((*nextNodeIter).first,  2); 
            bestEdgeConnect[refEnd] = refStart;
            bestEdgeConnect[altEnd] = altStart;
        }
    }
}

int VairiantGrpah::sumOfEdgePairQuality(PosAllele start, PosAllele end){
    std::map<int,VariantEdge*>::iterator edgeIter = edgeList.find( start.first );
    
    if( edgeIter != edgeList.end() ){
        if( start.second == 1 ){
            return edgeList[start.first]->ref->getQuality(end);
        }
        else if ( start.second == 2 ){
            return edgeList[start.first]->alt->getQuality(end);
        }
    }
    return -1;
}

VairiantGrpah::VairiantGrpah(std::string &in_ref, PhasingParameters &in_params){
    params=&in_params;
    ref=&in_ref;
}

VairiantGrpah::~VairiantGrpah(){
}

void VairiantGrpah::addEdge(std::vector<ReadVariant> &readVariant){
    // iter all read
    for(std::vector<ReadVariant>::iterator readIter = readVariant.begin() ; readIter != readVariant.end() ; readIter++ ){
        // record all variant 
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
            // this allele support ref
            if( (*variant1Iter).allele == 0 )
                edgeList[(*variant1Iter).position]->ref->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).read_name);
            // this allele support alt
            if( (*variant1Iter).allele == 1 )
                edgeList[(*variant1Iter).position]->alt->addSubEdge((*variant1Iter).quality, (*variant2Iter),(*readIter).read_name);
            variant1Iter++;
            variant2Iter++;
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

        for(std::map<int,VariantEdge*>::iterator nodeIter = edgeList.begin() ; nodeIter != edgeList.end() ; nodeIter++ ){
            std::vector<std::string> refEdge = (*nodeIter).second->ref->showEdge( "ref_" + std::to_string((*nodeIter).first) );
            std::vector<std::string> altEdge = (*nodeIter).second->alt->showEdge( "alt_" + std::to_string((*nodeIter).first) );
            
            for(std::vector<std::string>::iterator iter = refEdge.begin() ; iter != refEdge.end() ; iter++ )
                resultVcf << (*iter) << "\n";
            for(std::vector<std::string>::iterator iter = altEdge.begin() ; iter != altEdge.end() ; iter++ )
                resultVcf << (*iter) << "\n";
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

    // disjoint path process, find best edge
    this->findBestEdgeList();
    // find a path using current edge
    this->findResultPath();
    
    for(int round = 1 ; round <= 2 ; round++ ){
        // round 1: try to connect block by average quality
        // round 2: if the breakpoint with equal weighted. Use number of read connect two block.
        this->connectBlockByReadName(1, params->qualityThreshold, false, (round == 2) );
        this->findResultPath();

        // connect block by total quality
        this->connectBlockByQuality();
        this->findResultPath();
       
        // After the above steps, it is the result of the first phasing. 
        // Preliminary HP and Block results have been obtained. 
        // But the reason why many blocks are disconnected is because the choice of two weights is the same.
        // This step will directly use Read Name to determine whether the two blocks can be connected.
        // At this stage, a new edge will be added between two connectable blocks.
        for(int i = 1 ; i <= params->crossBlock ; i++ ){
            while(true){
                int beforeBlockSize = blockVec.size();
                bool connect = this->connectBlockByReadName(i, params->blockReadThreshold, true, false);
                // find a path using current edge
                this->findResultPath();
                int aftertBlockSize = blockVec.size();
                
                if(!connect || beforeBlockSize == aftertBlockSize){
                    break;
                }
            }
        }
    }
    
}


  