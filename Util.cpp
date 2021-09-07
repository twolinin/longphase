#include "Util.h"

void ReadVariant::sort(){
    std::sort(variantVec.begin(), variantVec.end(), less_than_key());
}

std::string getTargetString(std::string line, std::string start_sign, std::string end_sign){
    int start = line.find(start_sign) + 1;
    int end   = line.find(end_sign);
    int target_length = end - start;
    return line.substr(start,target_length);
}

int homopolimerLength(int snp_pos, const std::string &ref_string){
    int homopolimer_length = 1;
    char element = ref_string.at(snp_pos);
    int ref_len = ref_string.length();

    int pos = snp_pos-1;
    while( ref_string.at(pos) == element ){
        pos--;
        homopolimer_length++;
        if(homopolimer_length>=10)
            break;
    }

    pos = snp_pos+1;
    
    if( pos < ref_len ){
        while( ref_string.at(pos) == element ){
            pos++;
            homopolimer_length++;
            if( pos >= ref_len)
                break;
            if(homopolimer_length>=10)
            break;
        }
    }

    return homopolimer_length;
}







