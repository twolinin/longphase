make clean
make -j 32

#test_set1
  /bip8_disk/yungjen114/longphase/longphase phase \
  -s /big8_disk/giab_lsk114_2022.12/alignment/sup/pepper/pepper_sup_50x/OUTPUT_PREFIX.vcf.gz \
  -b /big8_disk/giab_lsk114_2022.12/alignment/sup/hg002.sup.50x.bam \
  -r /big8_disk/ref/GRCh38_no_alt_analysis_set.fasta \
  -t 32 \
  -o /big8_disk/yungjen114/longphase-output/longphase/quality_filter \
  --indels \
  --indelQuality 30 \
  --ont