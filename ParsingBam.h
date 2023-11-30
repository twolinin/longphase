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

#include <zlib.h>

struct RefAlt{
    std::string Ref;
    std::string Alt;
    bool is_reverse;
    bool is_modify;
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

class BaseVairantParser{
    public:
        // input parser
        void compressParser(std::string &variantFile);
        void unCompressParser(std::string &variantFile);
        virtual void parserProcess(std::string &input)=0;
        // output parser
        void compressInput(std::string variantFile, std::string resultFile, PhasingResult phasingResult);
        void unCompressInput(std::string variantFile, std::string resultFile, PhasingResult phasingResult);
        virtual void writeLine(std::string &input, bool &ps_def, std::ofstream &resultVcf, PhasingResult &phasingResult)=0;
};

class SnpParser : public BaseVairantParser{
    
    private:
        PhasingParameters *params;
        //std::string variantFile;
        // chr, variant position (0-base), allele haplotype
        std::map<std::string, std::map<int, RefAlt> > *chrVariant;
        // id and idx
        std::vector<std::string> chrName;
        // chr, variant position (0-base)
        std::map<std::string, std::map<int, bool> > chrVariantHomopolymer;
        
        // override input parser
        void parserProcess(std::string &input);
        // override output parser
        void writeLine(std::string &input, bool &ps_def, std::ofstream &resultVcf, PhasingResult &phasingResult);
        
        bool commandLine;
        
    public:
    
        SnpParser(PhasingParameters &in_params);
        ~SnpParser();
            
        std::map<int, RefAlt> getVariants(std::string chrName);  

        std::vector<std::string> getChrVec();
        
        bool findChromosome(std::string chrName);
        
        int getLastSNP(std::string chrName);
        
        void writeResult(PhasingResult phasingResult);

        bool findSNP(std::string chr, int posistion);
        
        void filterSNP(std::string chr, std::vector<ReadVariant> &readVariantVec, std::string &chr_reference);
};

class SVParser : public BaseVairantParser{
    
    private:
        PhasingParameters *params;
        SnpParser *snpFile;

        // chr , variant position (0-base), read
        std::map<std::string, std::map<int, std::map<std::string ,bool> > > *chrVariant;
        // chr, variant position (0-base)
        std::map<std::string, std::map<int, bool> > posDuplicate;
        
        // override input parser
        void parserProcess(std::string &input);
        // override output parser
        void writeLine(std::string &input, bool &ps_def, std::ofstream &resultVcf, PhasingResult &phasingResult);
        
        bool commandLine;
        
    public:
    
        SVParser(PhasingParameters &params, SnpParser &snpFile);
        ~SVParser();
            
        std::map<int, std::map<std::string ,bool> > getVariants(std::string chrName);  

        void writeResult(PhasingResult phasingResult);
};

class METHParser : public BaseVairantParser{
    
    private:
        PhasingParameters *params;
        SnpParser *snpFile;
        
        int representativePos;
        int upMethPos;
        
        // chr , variant position (0-base), read, (is reverse strand)
        std::map<std::string, std::map<int, std::map<std::string ,RefAlt> > > *chrVariant;
        // In a series of consecutive methylation positions, 
        // the first methylation position will be used as the representative after merging.
        // This map is used to find the coordinates that originally represented itself
        std::map<int, int > *representativeMap;

        // override input parser
        void parserProcess(std::string &input);
        // override output parser
        void writeLine(std::string &input, bool &ps_def, std::ofstream &resultVcf, PhasingResult &phasingResult);
        
        bool commandLine;
        
    public:
        
        std::map<int, std::map<std::string ,RefAlt> > getVariants(std::string chrName);  
        
        METHParser(PhasingParameters &params, SnpParser &snpFile);
        ~METHParser();
		
		void writeResult(PhasingResult phasingResult);
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
        std::vector<std::string> BamFileVec;
        // SNP map and iter
        std::map<int, RefAlt> *currentVariants;
        std::map<int, RefAlt>::iterator firstVariantIter;
        // SV map and iter
        std::map<int, std::map<std::string ,bool> > *currentSV;
        std::map<int, std::map<std::string ,bool> >::iterator firstSVIter;
        // mod map and iter
        std::map<int, std::map<std::string ,RefAlt> > *currentMod;
        std::map<int, std::map<std::string ,RefAlt> >::iterator firstModIter;
        void get_snp(const  bam_hdr_t &bamHdr,const bam1_t &aln, std::vector<ReadVariant> &readVariantVec, const std::string &ref_string, bool isONT);
   
    public:
        BamParser(std::string chrName, std::vector<std::string> inputBamFileVec, SnpParser &snpMap, SVParser &svFile, METHParser &modFile);
        ~BamParser();
        
        void direct_detect_alleles(int lastSNPPos, PhasingParameters params, std::vector<ReadVariant> &readVariantVec , const std::string &ref_string);

};


#endif