#include "Util.h"

void ReadVariant::sort(){
    std::sort(variantVec.begin(), variantVec.end(), less_than_key());
}

void mergeAllChrPhasingResult(const ChrPhasingResult& allChrPhasingResults, PhasingResult& mergedPhasingResult) {
    for(const auto& chrPair : allChrPhasingResults){
        const PhasingResult& singlePhasingResult = chrPair.second;
        mergedPhasingResult.insert(singlePhasingResult.begin(), singlePhasingResult.end());
    }
}

void setNumThreads(const int& defaultChrThreads,const int& availableThreads,  int& chrNumThreads, int& bamParserNumThreads){
    chrNumThreads = defaultChrThreads;
    bamParserNumThreads = 5;
    if(availableThreads > 0){
        if(availableThreads <= 24){
            chrNumThreads = availableThreads;
            bamParserNumThreads = 1;

        }else if(availableThreads <= 48){
            chrNumThreads = 12;
            bamParserNumThreads = std::max(1, availableThreads / 12);
        }else{
            chrNumThreads = availableThreads / 4;
            bamParserNumThreads = 4;
        }
    }
}

std::string getTargetString(std::string line, std::string start_sign, std::string end_sign){
    int start = line.find(start_sign) + 1;
    int end   = line.find(end_sign);
    int target_length = end - start;
    return line.substr(start,target_length);
}

int homopolymerLength(int snp_pos, const std::string &ref_string){
    int homopolymer_length = 1;
    int ref_len = ref_string.length();
    
    if( snp_pos + 1 >= ref_len ){
        return homopolymer_length;
    }
    
    char element = ref_string.at(snp_pos);
    
    int pos = snp_pos-1;
    while( ref_string.at(pos) == element ){
        pos--;
        homopolymer_length++;
        if(homopolymer_length>=10 || pos < 0 )
            break;
    }

    pos = snp_pos+1;
    
    if( pos < ref_len ){
        while( ref_string.at(pos) == element ){
            pos++;
            homopolymer_length++;
            if( pos >= ref_len)
                break;
            if(homopolymer_length>=10)
            break;
        }
    }

    return homopolymer_length;
}