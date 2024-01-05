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

void setPhasingNumThreads(const int& defaultChrThreads,const int& availableThreads,  int& chrNumThreads, int& bamParserNumThreads){
    if(availableThreads == 0){
        // The default thread value is 0, indicating that each chromosome concurrently utilizes 
        // a single thread for computation. Due to the limited multithreading acceleration space 
        // in htslib and an average utilization of no more than 5 threads per chromosome, the 
        // default value is set to 5.
        chrNumThreads = defaultChrThreads;
        bamParserNumThreads = 5;
    }
    else if( defaultChrThreads > availableThreads){
        // Due to the number of available threads being less than the number of chromosomes, 
        // processing is performed on a subset of chromosomes equal to the availableThreads count, 
        // with a limitation on the number of threads for reading the BAM file for each chromosome.
        chrNumThreads = availableThreads;
        bamParserNumThreads = 1;
    }
    else{
        // If the number of available threads exceeds the count of chromosomes, the surplus 
        // threads will be utilized to accelerate the speed of htslib in parsing BAM files.
        chrNumThreads = defaultChrThreads;
        bamParserNumThreads = std::max(1,availableThreads/defaultChrThreads);
    }
}

void setModcallNumThreads(const int& availableThreads,  int& chrNumThreads, int& bamParserNumThreads){
    // Using 4 threads in htslib tends to yield higher utilization.
    // Allocate threads to the BAM parser, but do not exceed the maximum number of threads.
    // This ensures that, even when more threads are available, the BAM parser does not use too many,
    // thereby leaving resources available for other tasks.
    const int maxBamParserThreads = 4;
    bamParserNumThreads = std::min(maxBamParserThreads, availableThreads);

    // Allocate at least one thread to the chromosome processing task.
    // If there are enough available threads, the remaining threads will be allocated to this task.
    // This allocation maintains a balance in the total number of threads used.
    chrNumThreads = std::max(1, availableThreads / bamParserNumThreads);
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