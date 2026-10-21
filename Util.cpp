#include "Util.h"
#include <condition_variable>
#include <mutex>
#include <thread>
#include <htslib/kseq.h>

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
    // pos is -1 for a variant at the first base
    while( pos >= 0 && ref_string.at(pos) == element ){
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

namespace {

// Threads that run a job, fn(t) for t = 0 .. n-1 (0 on the calling thread),
// each time run() is called, and sleep in between (OpenMP threads would spin
// while the batches are read and written).
class BatchThreads {
public:
    explicit BatchThreads(int n) : n_(n) {
        for (int t = 1; t < n; ++t)
            threads_.emplace_back([this, t]() {
                uint64_t seen = 0;
                while (true) {
                    const std::function<void(int)>* job;
                    {
                        std::unique_lock<std::mutex> lock(m_);
                        cv_.wait(lock, [&]() { return stop_ || gen_ != seen; });
                        if (stop_) return;
                        seen = gen_;
                        job = job_;
                    }
                    (*job)(t);
                    std::lock_guard<std::mutex> lock(m_);
                    if (--pending_ == 0) doneCv_.notify_one();
                }
            });
    }
    ~BatchThreads() {
        {
            std::lock_guard<std::mutex> lock(m_);
            stop_ = true;
        }
        cv_.notify_all();
        for (auto& th : threads_) th.join();
    }
    void run(const std::function<void(int)>& fn) {
        {
            std::lock_guard<std::mutex> lock(m_);
            job_ = &fn;
            pending_ = n_ - 1;
            ++gen_;
        }
        cv_.notify_all();
        fn(0);
        std::unique_lock<std::mutex> lock(m_);
        doneCv_.wait(lock, [&]() { return pending_ == 0; });
    }
private:
    int n_;
    std::vector<std::thread> threads_;
    std::mutex m_;
    std::condition_variable cv_, doneCv_;
    const std::function<void(int)>* job_ = nullptr;
    uint64_t gen_ = 0;
    int pending_ = 0;
    bool stop_ = false;
};

}

void readVcfRecords(htsFile* fp, bcf_hdr_t* hdr, int threads, size_t batch,
                    const VcfRecordFn& work, const VcfRecordFn& done) {
    threads = std::max(1, threads);
    std::vector<bcf_hdr_t*> hdrs(threads, hdr);
    for (int t = 1; t < threads; ++t) hdrs[t] = bcf_hdr_dup(hdr);
    std::vector<kstring_t> lines(batch);
    for (auto& l : lines) l = KS_INITIALIZE;
    std::vector<bcf1_t*> recs(batch);
    for (auto& r : recs) r = bcf_init();
    std::vector<int> parsed(batch), owner(batch);
    BatchThreads pool(threads);
    size_t n = 0;
    // thread t parses records [n * t / threads, n * (t + 1) / threads)
    const std::function<void(int)> job = [&](int t) {
        for (size_t i = n * t / threads; i < n * (t + 1) / threads; ++i) {
            owner[i] = t;
            parsed[i] = vcf_parse(&lines[i], hdrs[t], recs[i]);
            if (parsed[i] == 0) work(i, hdrs[t], recs[i]);
        }
    };
    bool more = true;
    while (more) {
        n = 0;
        while (n < batch && hts_getline(fp, KS_SEP_LINE, &lines[n]) >= 0) ++n;
        if (n < batch) more = false;
        if (n == 0) break;
        pool.run(job);
        for (size_t i = 0; i < n; ++i) {
            if (parsed[i] != 0) { more = false; break; }
            done(i, hdrs[owner[i]], recs[i]);
        }
    }
    for (auto& r : recs) bcf_destroy(r);
    for (auto& l : lines) free(l.s);
    for (int t = 1; t < threads; ++t) bcf_hdr_destroy(hdrs[t]);
}
