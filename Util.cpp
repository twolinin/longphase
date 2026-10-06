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

std::string getTargetString(std::string line, std::string start_sign, std::string end_sign){
    int start = line.find(start_sign) + 1;
    int end   = line.find(end_sign);
    int target_length = end - start;
    return line.substr(start,target_length);
}

int findFieldValueStart(const std::string &format, const std::string &sample,
                        const std::string &tag) {
  int colon_pos = 0;
  int tag_pos = format.find(tag);
  for (int i = 0; i < tag_pos; i++) {
    if (format[i] == ':')
      colon_pos++;
  }
  int current_colon = 0;
  int value_start = 0;
  for (unsigned int i = 0; i < sample.length(); i++) {
    if (current_colon >= colon_pos)
      break;
    if (sample[i] == ':')
      current_colon++;
    value_start++;
  }
  return value_start;
}

bool isDangerIndel(int variant_pos, const std::string &ref, size_t ref_allele_len,
                   size_t alt_allele_len) {
  if (ref_allele_len <= 1 && alt_allele_len <= 1)
    return false;

  std::string repeat = ref.substr(variant_pos + 1, 2);
  int ref_pos = variant_pos;
  for (int i = 0; i < 5; i++) {
    if (repeat[0] != ref[ref_pos + 1] || repeat[1] != ref[ref_pos + 2])
      return false;
    ref_pos += 2;
  }
  return true;
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

template <typename T>
double calculateMean(const std::vector<T>& data) {
    if(data.empty()) {
        return 0;
    }
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    return sum / data.size();
}
template double calculateMean<int>(const std::vector<int>& data);
template double calculateMean<double>(const std::vector<double>& data);