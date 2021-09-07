#ifndef PARSINGBAM_H
#define PARSINGBAM_H

#include "Util.h"
#include "PhasingProcess.h"
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/kbitset.h>
#include <htslib/thread_pool.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

struct RefAlt{
    std::string Ref;
    std::string Alt;
};

class FastaParser{
    private:
        // file name
        std::string fastaFile ;
        std::vector<std::string> chrName;
        std::vector<int> last_pos;
    public:
        FastaParser(std::string fastaFile  , std::vector<std::string> chrName , std::vector<int> last_pos);
        ~FastaParser();
        
        // chrName, chr string
        std::map<std::string, std::string > chrString;
    
};

class VariantParser{
    
    private:
    
        std::string variantFile;
        // chr, variant position (0-base), allele haplotype
        std::map<std::string, std::map<int, RefAlt> > chrVariant;
        // id and idx
        std::vector<std::string> chrName;
        // chr, variant position (0-base)
        std::map<std::string, std::map<int, bool> > chrVariantHomopolymer;
        
        int homopolimerEnd(int snp_pos, const std::string &ref_string);
    
    public:
    
        VariantParser(std::string variantFile);
        ~VariantParser();
            
        std::map<int, RefAlt> getVariants(std::string chrName);  
        
        std::map<int, bool> getHomopolymeVariants(std::string chrName);  
        
        std::vector<std::string> getChrVec();
        
        bool findChromosome(std::string chrName);
        
        int getLastSNP(std::string chrName);
        
        void writeResult(std::string resultPrefix, PhasingResult phasingResult);
        
        bool findSNP(std::string chr, int posistion);
        
        void filterSNP(std::string chr, std::vector<ReadVariant> &readVariantVec);
};

class SVParser{
    
    private:
    
        std::string variantFile;
        // chr , variant position (0-base), read
        std::map<std::string, std::map<int, std::map<std::string ,bool> > > chrVariant;
        // chr, variant position (0-base)
        std::map<std::string, std::map<int, bool> > posDuplicate;

    public:
    
        SVParser(std::string variantFile, VariantParser snpFile);
        ~SVParser();
            
        std::map<int, std::map<std::string ,bool> > getVariants(std::string chrName);  

        void writeResult(std::string resultPrefix, PhasingResult phasingResult);
};

struct Alignment{
    std::string chr;
    std::string qname;
    int refStart;
    int qlen;
    char *qseq;
    int cigar_len;
    int *op;
    int *ol;
    char *quality;
    bool is_reverse;
};

class BamParser{
    
    private:
        std::string chrName;
        std::string BamFile;
        // SNP map and iter
        std::map<int, RefAlt> currentVariants;
        std::map<int, RefAlt>::iterator firstVariantIter;
        // SV map and iter
        std::map<int, std::map<std::string ,bool> > currentSV;
        std::map<int, std::map<std::string ,bool> >::iterator firstSVIter;
        void get_snp(Alignment align, std::vector<ReadVariant> &readVariantVec, const std::string &ref_string, bool isONT);
        //int homopolimerLength(int snp_pos, const std::string &ref_string);
        bool continueHomopolimer(int snp_pos, const std::string &ref_string);
    
    public:
        BamParser(std::string chrName, std::string inputBamFile, VariantParser snpMap, SVParser svFile);
        ~BamParser();
        
        void direct_detect_alleles(int lastSNPPos, PhasingParameters params, std::vector<ReadVariant> &readVariantVec , const std::string &ref_string);

};


#endif